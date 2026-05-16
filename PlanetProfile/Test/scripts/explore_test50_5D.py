"""5D variant of Test50 MCMC — thin shim that delegates to Test50 with the 5D config.

Drops Tb_K from sampling (fixed at 250.965 K, the upper grid edge) and collapses
HP-ice viscosities (Ice III/V/VI) into a single log10_eta_HP via param_groups.

Param space: [alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil]

To run:
  mamba activate PPcl
  python PlanetProfile/Test/scripts/explore_test50_5D.py
"""
import os
import sys
import runpy

_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_config_5D = os.path.join(_pp_root, 'PlanetProfile', 'Inference', 'configs',
                           'test50_titan_noocean_andrade_5D.json')

# Delegate to Test50 main with the 5D config
sys.argv = [sys.argv[0], '--config', _config_5D]
runpy.run_path(
    os.path.join(_pp_root, 'PlanetProfile', 'Test',
                 'Test50_mcmc_andrade_noocean_yao2014.py'),
    run_name='__main__'
)
