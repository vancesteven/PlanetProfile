import sys
import os
import json
import numpy as np
import logging
from typing import Dict, Any, Tuple

# Path setup
_pp_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _pp_root)
# Add era implementation to path
era_impl_path = os.path.join(_pp_root, 'era', 'implementation')
if os.path.exists(era_impl_path):
    sys.path.insert(0, era_impl_path)

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

# Attempt to import era; fallback if not available
try:
    from futs import Problem, Solution, search
    HAS_ERA = True
except ImportError:
    HAS_ERA = False
    print("Warning: 'era' tool implementation not found in path.")

# --- PlanetProfile Era Integration ---

def execute_titan_config(problem: Any, solution: Any) -> float:
    """
    Execute function for 'era' search.
    Parses configuration, runs MCMC, and scores the result.
    """
    try:
        config_data = json.loads(solution.program)
        config = InferenceConfig.from_dict(config_data)
        
        # Ensure we run in a mode that's fast but informative
        config.sampler_settings['n_effective'] = 200
        config.sampler_settings['n_reeval'] = 50
        
        runner = MCMCRunner(config)
        result = runner.run()
        
        # Scoring Logic:
        # 1. Convergence (ESS) - max 10 pts
        ess = result.convergence_metrics.get('ess', 0)
        ess_score = min(10.0, ess / 20.0)
        
        # 2. Data Fit (max log-likelihood) - max 10 pts
        # Best possible chi2=0 -> logl=0. Penalty for large chi2.
        max_logl = np.max(result.log_likelihoods)
        fit_score = max(0.0, 10.0 + max_logl)
        
        # 3. Physical Plausibility (Heat Partitioning) - max 10 pts
        # For Titan, we strongly expect dissipation in the ice layers (Ih or HP).
        # Solutions where 100% of heat is in the silicate are likely unphysical 
        # or indicate the sampler is stuck in a bad local minimum.
        heating = result.heating_results
        if heating:
            sil_fracs = [h.get('Silicate_W', 0) / (sum(v for k, v in h.items() if not k.endswith('_W')) + 1e-30) 
                         for h in heating]
            med_sil_frac = np.median(sil_fracs)
            # Higher score for lower silicate fraction (preferring ice dissipation)
            heat_score = 10.0 * (1.0 - med_sil_frac)
        else:
            heat_score = 0.0
            
        total_score = ess_score + fit_score + heat_score
        print(f"  [Era] Eval complete. Score: {total_score:.2f} (ESS: {ess_score:.1f}, Fit: {fit_score:.1f}, Heat: {heat_score:.1f})")
        return total_score
        
    except Exception as e:
        print(f"  [Era] Evaluation failed: {e}")
        return 0.0

def generate_improved_config(problem: Any, parent_solution: Any, parent_score: float) -> Any:
    """
    Generator function for 'era' search.
    In a production setup, this would call the LLM with a detailed prompt.
    """
    # For now, this is a placeholder that would be wired to the LLM.
    # The LLM would be prompted with:
    # "The current configuration has a score of {parent_score}. 
    #  JSON: {parent_solution.program}
    #  Improve the prior bounds to better capture the posterior and improve convergence."
    return parent_solution

def run_era_optimization(num_iterations: int = 5):
    """Main entry point for Titan MCMC meta-optimization."""
    if not HAS_ERA:
        print("Error: 'era' tool is required for this optimization.")
        return

    # Load initial solution (standard Test50 8D)
    config_path = os.path.join(_pp_root, 'PlanetProfile', 'Inference', 'configs', 'test50_titan_noocean_andrade_8D.json')
    with open(config_path, 'r') as f:
        initial_json = f.read()
    
    problem = Problem(description="Find optimal MCMC priors for No-Ocean Titan to maximize convergence and physical plausibility.")
    initial_sol = Solution(program=initial_json)
    
    print(f"Starting Era FUTS search for {num_iterations} iterations...")
    best_sol, best_score = search(
        problem=problem,
        initial_solution=initial_sol,
        initial_score=0.0, # Will be re-evaluated
        generate_fn=generate_improved_config,
        execute_fn=execute_titan_config,
        num_iterations=num_iterations
    )
    
    print("\n" + "="*50)
    print(f"Best Score Found: {best_score:.2f}")
    print("Best Configuration:")
    print(best_sol.program)
    print("="*50)

if __name__ == "__main__":
    # Note: Requires a functional PP environment and 'pocomc'
    # run_era_optimization()
    print("Titan Era Optimization Script ready.")
    print("To run, uncomment run_era_optimization() in __main__.")
