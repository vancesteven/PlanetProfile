"""
Diagnostic script to test MCMC forward model at prior center.

Tests:
1. Load structure cache
2. Call forward_model_k2 with default Andrade parameters
3. Check if k2 is NaN or valid
4. Compute log-likelihood value

Author: Diagnostic script
Date: 2026-04-29
"""
import numpy as np
import sys
from pathlib import Path

# Add PlanetProfile to path
sys.path.insert(0, str(Path(__file__).parent))

from PlanetProfile.Inference.structure_cache import load_structure_cache
from PlanetProfile.Inference.forward_models import forward_model_k2, create_log_likelihood
from PlanetProfile.Inference.parameter_registry import PARAMETER_PRESETS, PARAMETER_REGISTRY

# ============================================================================
# Test Configuration
# ============================================================================

# Use andrade_titan preset
preset = PARAMETER_PRESETS['andrade_titan']
param_ids = preset['parameters']
observables = preset['observables']

print("=" * 70)
print("MCMC FORWARD MODEL DIAGNOSTIC")
print("=" * 70)
print(f"\nPreset: {preset['name']}")
print(f"Parameters: {param_ids}")
print(f"Observables: {observables}")

# Get default parameter values (prior centers)
theta_default = []
theta_dict = {}
for param_id in param_ids:
    param_def = PARAMETER_REGISTRY[param_id]

    if param_def.default_prior in ['uniform', 'log-uniform']:
        # Use midpoint of bounds
        low, high = param_def.default_bounds
        center = (low + high) / 2.0
    elif param_def.default_prior == 'normal':
        center = param_def.default_mean
    else:
        center = 0.0

    theta_default.append(center)
    theta_dict[param_id] = center
    print(f"  {param_id} = {center:.4f} ({param_def.label})")

theta_default = np.array(theta_default)

# ============================================================================
# Load Structure Cache
# ============================================================================

# Try to find Titan structure cache
cache_candidates = [
    'titan_structure.pkl',
    'titan_cache/titan_structure_noclath.pkl',
    'titan_cache/titan_structure_clath.pkl',
    'Titan/structures_titan_noclath.pkl',
    'Test41/structures_titan_noclath.pkl',
    'structures_titan_noclath.pkl',
]

structure_data = None
cache_path = None

for candidate in cache_candidates:
    candidate_path = Path(candidate)
    if candidate_path.exists():
        cache_path = candidate_path
        break

    # Try in current directory
    candidate_path = Path(__file__).parent / candidate
    if candidate_path.exists():
        cache_path = candidate_path
        break

if cache_path is None:
    print("\n❌ ERROR: Could not find structure cache file!")
    print("   Tried:")
    for c in cache_candidates:
        print(f"     - {c}")
    print("\n   Generate cache with: python -m PlanetProfile.Inference.prepare_structure_variants Test41")
    sys.exit(1)

print(f"\n✅ Loading structure cache: {cache_path}")
structure_data = load_structure_cache(str(cache_path), validate_bodyname=None)

print(f"   Loaded {structure_data['n_layers']} layers")
print(f"   Phases: {structure_data['region_phases']}")
print(f"   Body: {structure_data.get('bodyname', 'unknown')}")

# ============================================================================
# Test 1: Single Forward Model Call
# ============================================================================

print("\n" + "=" * 70)
print("TEST 1: Forward Model Call at Prior Center")
print("=" * 70)

Re_k2, Im_k2, _ = forward_model_k2(
    theta_default,
    structure_data,
    rheology='andrade',
    return_heating=False,
    arrhenius_params=None
)

print(f"\nResults:")
print(f"  Re(k₂) = {Re_k2:.6f}")
print(f"  Im(k₂) = {Im_k2:.6f}")
print(f"  |k₂| = {np.sqrt(Re_k2**2 + Im_k2**2):.6f}")

if np.isnan(Re_k2):
    print("\n❌ ERROR: Forward model returned NaN!")
    print("   TidalPy solver likely failed.")
    sys.exit(1)
else:
    print("\n✅ Forward model succeeded (no NaN)")

# ============================================================================
# Test 2: Log-Likelihood at Prior Center
# ============================================================================

print("\n" + "=" * 70)
print("TEST 2: Log-Likelihood at Prior Center")
print("=" * 70)

log_likelihood_fn = create_log_likelihood(
    observables,
    structure_data,
    rheology='andrade',
    arrhenius_params=None
)

log_like_center = log_likelihood_fn(theta_default)

print(f"\nLog-likelihood: {log_like_center:.6f}")

if log_like_center == -1e30:
    print("❌ Log-likelihood = -1e30 (TidalPy failed)")
    sys.exit(1)
elif log_like_center < -1e10:
    print("⚠️  WARNING: Extremely low log-likelihood!")
    print("   This indicates very poor fit at prior center.")
    print("   MCMC may struggle to find high-likelihood regions.")
elif np.isinf(log_like_center):
    print("❌ Log-likelihood is -inf!")
    sys.exit(1)
else:
    print("✅ Log-likelihood is finite")

# Compute chi-squared
chi2 = -2.0 * log_like_center
print(f"χ² = {chi2:.2f}")

# Check against observations
print(f"\nObservations vs. Model:")
for obs_name, (obs_val, obs_err) in observables.items():
    if obs_name == 'Re_k2':
        model_val = Re_k2
    elif obs_name == 'Im_k2':
        model_val = Im_k2
    else:
        continue

    residual = model_val - obs_val
    sigma_deviation = residual / obs_err

    print(f"  {obs_name}:")
    print(f"    Observed: {obs_val:.4f} ± {obs_err:.4f}")
    print(f"    Model:    {model_val:.4f}")
    print(f"    Residual: {residual:.4f} ({sigma_deviation:.2f}σ)")

# ============================================================================
# Test 3: Sample Random Points from Prior
# ============================================================================

print("\n" + "=" * 70)
print("TEST 3: Sample 10 Random Points from Prior")
print("=" * 70)

from scipy.stats import uniform, loguniform

n_test = 10
success_count = 0
nan_count = 0
low_like_count = 0

for i in range(n_test):
    # Sample from prior
    theta_sample = []
    for param_id in param_ids:
        param_def = PARAMETER_REGISTRY[param_id]

        if param_def.default_prior == 'uniform':
            low, high = param_def.default_bounds
            val = np.random.uniform(low, high)
        elif param_def.default_prior == 'log-uniform':
            low, high = param_def.default_bounds
            val = loguniform(a=10**low, b=10**high).rvs()
            val = np.log10(val)  # Convert back to log-space
        elif param_def.default_prior == 'normal':
            val = np.random.normal(param_def.default_mean, param_def.default_std)
        else:
            val = 0.0

        theta_sample.append(val)

    theta_sample = np.array(theta_sample)

    # Evaluate
    log_like = log_likelihood_fn(theta_sample)

    if log_like == -1e30:
        nan_count += 1
        status = "NaN"
    elif log_like < -100:
        low_like_count += 1
        status = "Low"
    else:
        success_count += 1
        status = "OK"

    print(f"  Sample {i+1}: log(L) = {log_like:10.2f}  [{status}]")

print(f"\nSummary:")
print(f"  Success: {success_count}/{n_test}")
print(f"  NaN:     {nan_count}/{n_test}")
print(f"  Low:     {low_like_count}/{n_test}")

if success_count == 0:
    print("\n❌ CRITICAL: All prior samples failed or have low likelihood!")
    print("   MCMC cannot explore if forward model consistently fails.")
elif success_count < 5:
    print("\n⚠️  WARNING: Most prior samples have low/failed likelihood.")
    print("   MCMC may have very low acceptance rate.")
else:
    print("\n✅ Sufficient prior samples succeeded.")

print("\n" + "=" * 70)
print("DIAGNOSTIC COMPLETE")
print("=" * 70)
