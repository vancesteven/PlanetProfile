"""
Diagnose acceptance rate issue by evaluating log-likelihood at known good parameters.

The PPTest41 sanity check uses theta = [0.2, log10(0.05), 12.0, 13.0, 19.0]
which produces k2 ≈ 0.280 - 0.144j. This should give reasonable likelihood
against observations Re_k2=0.608±0.048, abs_Im_k2=0.135±0.035.

Author: Diagnostic script
Date: 2026-04-29
"""
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from PlanetProfile.Inference.structure_cache import load_structure_cache
from PlanetProfile.Inference.forward_models import forward_model_k2, create_log_likelihood

print("=" * 70)
print("ACCEPTANCE RATE DIAGNOSTIC")
print("=" * 70)

# Load structure cache
cache_path = Path('titan_structure.pkl')
if not cache_path.exists():
    print(f"\n❌ ERROR: {cache_path} not found")
    sys.exit(1)

print(f"\n✅ Loading: {cache_path}")
structure_data = load_structure_cache(str(cache_path), validate_bodyname=None)

# Define observables (as used in presets)
observables = {
    'Re_k2': (0.608, 0.048),
    'abs_Im_k2': (0.135, 0.035)
}

print("\nObservables:")
print(f"  Re(k₂) = {observables['Re_k2'][0]:.3f} ± {observables['Re_k2'][1]:.3f}")
print(f"  |Im(k₂)| = {observables['abs_Im_k2'][0]:.3f} ± {observables['abs_Im_k2'][1]:.3f}")

# Create log-likelihood function
log_likelihood_fn = create_log_likelihood(
    observables,
    structure_data,
    rheology='andrade',
    arrhenius_params=None
)

# ============================================================================
# Test 1: PPTest41 Sanity Check Parameters
# ============================================================================

print("\n" + "=" * 70)
print("TEST 1: PPTest41 Sanity Check Parameters")
print("=" * 70)

theta_sanity = np.array([
    0.2,              # alpha
    np.log10(0.05),   # log10_zeta ≈ -1.3
    12.0,             # log10_eta_Ih
    13.0,             # log10_eta_HP
    19.0              # log10_eta_sil
])

print("\nParameters:")
print(f"  alpha = {theta_sanity[0]:.3f}")
print(f"  log10(zeta) = {theta_sanity[1]:.3f}  (zeta = {10**theta_sanity[1]:.4f})")
print(f"  log10(eta_Ih) = {theta_sanity[2]:.1f}  (eta_Ih = {10**theta_sanity[2]:.1e} Pa·s)")
print(f"  log10(eta_HP) = {theta_sanity[3]:.1f}  (eta_HP = {10**theta_sanity[3]:.1e} Pa·s)")
print(f"  log10(eta_sil) = {theta_sanity[4]:.1f} (eta_sil = {10**theta_sanity[4]:.1e} Pa·s)")

# Forward model
Re_k2, Im_k2, _ = forward_model_k2(
    theta_sanity,
    structure_data,
    rheology='andrade',
    return_heating=False,
    arrhenius_params=None
)

print(f"\nForward Model Output:")
print(f"  Re(k₂) = {Re_k2:.6f}")
print(f"  Im(k₂) = {Im_k2:.6f}")
print(f"  |Im(k₂)| = {abs(Im_k2):.6f}")
print(f"  |k₂| = {np.sqrt(Re_k2**2 + Im_k2**2):.6f}")

# Log-likelihood
log_like = log_likelihood_fn(theta_sanity)
chi2 = -2.0 * log_like

print(f"\nLog-Likelihood: {log_like:.6f}")
print(f"χ² = {chi2:.2f}")

# Residuals
Re_residual = Re_k2 - observables['Re_k2'][0]
Re_sigma = Re_residual / observables['Re_k2'][1]

Im_residual = abs(Im_k2) - observables['abs_Im_k2'][0]
Im_sigma = Im_residual / observables['abs_Im_k2'][1]

print(f"\nResiduals:")
print(f"  Re(k₂):  {Re_residual:+.4f}  ({Re_sigma:+.2f}σ)")
print(f"  |Im(k₂)|: {Im_residual:+.4f}  ({Im_sigma:+.2f}σ)")

# ============================================================================
# Test 2: Try Parameters Closer to Observations
# ============================================================================

print("\n" + "=" * 70)
print("TEST 2: Parameters Optimized Toward Observations")
print("=" * 70)

# Try parameters that might produce higher k2
# Lower viscosities -> more dissipation -> higher k2
theta_optimized = np.array([
    0.25,             # alpha (middle of prior)
    np.log10(0.1),    # log10_zeta ≈ -1.0 (higher zeta)
    13.5,             # log10_eta_Ih (middle)
    14.0,             # log10_eta_HP (middle)
    19.5              # log10_eta_sil (middle)
])

print("\nParameters (trial):")
print(f"  alpha = {theta_optimized[0]:.3f}")
print(f"  log10(zeta) = {theta_optimized[1]:.3f}  (zeta = {10**theta_optimized[1]:.4f})")
print(f"  log10(eta_Ih) = {theta_optimized[2]:.1f}")
print(f"  log10(eta_HP) = {theta_optimized[3]:.1f}")
print(f"  log10(eta_sil) = {theta_optimized[4]:.1f}")

Re_k2_opt, Im_k2_opt, _ = forward_model_k2(
    theta_optimized,
    structure_data,
    rheology='andrade',
    return_heating=False,
    arrhenius_params=None
)

log_like_opt = log_likelihood_fn(theta_optimized)
chi2_opt = -2.0 * log_like_opt

print(f"\nForward Model Output:")
print(f"  Re(k₂) = {Re_k2_opt:.6f}")
print(f"  Im(k₂) = {Im_k2_opt:.6f}")
print(f"  |Im(k₂)| = {abs(Im_k2_opt):.6f}")

print(f"\nLog-Likelihood: {log_like_opt:.6f}")
print(f"χ² = {chi2_opt:.2f}")

# ============================================================================
# Test 3: Grid Search to Find High-Likelihood Region
# ============================================================================

print("\n" + "=" * 70)
print("TEST 3: Grid Search for High-Likelihood Parameters")
print("=" * 70)

print("\nSearching parameter space (may take 30-60 seconds)...")

best_log_like = -np.inf
best_theta = None
best_k2 = None

n_samples = 0
high_like_count = 0  # Count samples with log_like > -10

# Grid over key parameters
for alpha in [0.2, 0.25, 0.3, 0.35]:
    for log10_zeta in [-1.5, -1.0, -0.5, 0.0]:
        for log10_eta_Ih in [12.0, 13.0, 14.0, 15.0]:
            theta_trial = np.array([alpha, log10_zeta, log10_eta_Ih, 14.0, 19.0])

            log_like_trial = log_likelihood_fn(theta_trial)
            n_samples += 1

            if log_like_trial > -10:
                high_like_count += 1

            if log_like_trial > best_log_like:
                best_log_like = log_like_trial
                best_theta = theta_trial.copy()

                # Get k2 for best
                Re_best, Im_best, _ = forward_model_k2(
                    best_theta, structure_data, rheology='andrade',
                    return_heating=False, arrhenius_params=None
                )
                best_k2 = (Re_best, Im_best)

print(f"\nGrid search complete: {n_samples} samples evaluated")
print(f"High-likelihood samples (log(L) > -10): {high_like_count}/{n_samples}")

print(f"\nBest parameters found:")
print(f"  alpha = {best_theta[0]:.3f}")
print(f"  log10(zeta) = {best_theta[1]:.3f}")
print(f"  log10(eta_Ih) = {best_theta[2]:.1f}")
print(f"  log10(eta_HP) = {best_theta[3]:.1f}")
print(f"  log10(eta_sil) = {best_theta[4]:.1f}")

print(f"\nBest k₂:")
print(f"  Re(k₂) = {best_k2[0]:.6f}")
print(f"  Im(k₂) = {best_k2[1]:.6f}")

print(f"\nBest log-likelihood: {best_log_like:.6f}")
print(f"Best χ² = {-2.0*best_log_like:.2f}")

# ============================================================================
# Analysis
# ============================================================================

print("\n" + "=" * 70)
print("DIAGNOSIS")
print("=" * 70)

if best_log_like < -50:
    print("\n❌ PROBLEM: Even best parameters have very low likelihood (log(L) < -50)")
    print("   This explains 0% acceptance rate — sampler is stuck in likelihood desert.")
    print("\n   Possible causes:")
    print("   1. Prior bounds don't cover high-k2 parameter region")
    print("   2. Single-zeta model fundamentally can't reach k2 ~ 0.6")
    print("   3. Observational uncertainties too small (over-constraining)")
elif best_log_like < -10:
    print("\n⚠️  Best log(L) ≈ {:.1f} indicates moderate fit possible".format(best_log_like))
    print("   MCMC might work with very long chains (low acceptance).")
else:
    print(f"\n✅ Found high-likelihood region (log(L) = {best_log_like:.1f})")
    print("   MCMC should work if starting near these parameters.")

print("\n" + "=" * 70)
