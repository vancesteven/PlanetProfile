
import pickle
import numpy as np
import os
import sys

# Environment setup
_pp_root = os.getcwd()
sys.path.insert(0, _pp_root)

from PlanetProfile.Test.Test49_mcmc_andrade_yao2014_clathrate import forward_model

def analyze_titan_sensitivity():
    # Load Test 49 (4km Clath)
    with open('PlanetProfile/Test/mcmc_results/Titan/Test49_clathrate4km/test49_clathrate4km_mcmc_results.pkl', 'rb') as f:
        res49 = pickle.load(f)

    # Load Test 48 (5km Clath / Path B)
    with open('PlanetProfile/Test/mcmc_results/Titan/Test48_andrade_yao2014/hybrid_hydro_andrade_yao2014_mcmc.pkl', 'rb') as f:
        res48 = pickle.load(f)

    # Load Cache 49
    with open('PlanetProfile/Test/mcmc_results/Titan/Test49_clathrate4km/titan_yao2014_clathrate4km_hybrid_hydro_grid_cache.pkl', 'rb') as f:
        cache49_data = pickle.load(f)
    grid_cache49 = cache49_data['grid_cache']
    tb_vals49 = np.array(sorted(set(k[0] for k in grid_cache49)))
    d_vals49 = np.array(sorted(set(k[1] for k in grid_cache49)))

    def get_stats(samples, params):
        return {p: (np.median(samples[:,i]), np.std(samples[:,i])) for i, p in enumerate(params)}

    stats49 = get_stats(res49['samples'], res49['params'])
    stats48 = get_stats(res48['samples'], res48['param_names'])

    print("--- PARAMETER COMPARISON (Median ± Std) ---")
    print(f"{'Param':<15} | {'Test 49 (4km)':<20} | {'Test 48 (5km)':<20}")
    print("-" * 60)
    for p in res49['params']:
        v49 = stats49[p]
        v48 = stats48.get(p, (np.nan, np.nan))
        print(f"{p:<15} | {v49[0]:>8.3f} ± {v49[1]:<8.3f} | {v48[0]:>8.3f} ± {v48[1]:<8.3f}")

    # Calculate heating for 49
    print("\n--- HEATING DISTRIBUTION (Test 49) ---")
    samples49 = res49['samples']
    n_eval = 100
    rng = np.random.default_rng(42)
    idx = rng.choice(len(samples49), n_eval, replace=False)
    heats = []
    for i in idx:
        _, _, _, _, h = forward_model(samples49[i], grid_cache49, tb_vals49, d_vals49, return_heating=True)
        if h: heats.append(h)
    
    if heats:
        totals = np.array([sum(h.values()) for h in heats])
        sil_frac = np.array([h.get('Sil', 0) for h in heats]) / totals
        ih_frac = np.array([h.get('Ih', 0) for h in heats]) / totals
        hp_frac = np.array([h.get('III', 0)+h.get('V', 0)+h.get('VI', 0) for h in heats]) / totals
        
        print(f"Median Silicate heating: {np.median(sil_frac):.1%}")
        print(f"Median Ice Ih heating:   {np.median(ih_frac):.1%}")
        print(f"Median HP Ice heating:   {np.median(hp_frac):.1%}")
    else:
        print("No heating data available.")

if __name__ == "__main__":
    analyze_titan_sensitivity()
