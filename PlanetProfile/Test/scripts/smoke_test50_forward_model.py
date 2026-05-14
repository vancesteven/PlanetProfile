"""Smoke test for Test50 8D forward_model.

Loads the rebuilt structure grid cache and exercises forward_model on a small
set of parameter vectors covering the prior corners + a no-ocean-safeguard
trigger case.  Used to validate the 8D pipeline (Tb sampling + grid interp +
no-ocean safeguard + Arrhenius Ih) end-to-end before launching the full
pocoMC sampler.
"""
import logging
import sys
import time

import numpy as np

logging.disable(logging.WARNING)

sys.path.insert(0, '.')
from PlanetProfile.Test.Test50_mcmc_andrade_noocean_yao2014 import (
    build_or_load_structure_grid, forward_model, PARAM_NAMES,
)


CASES = [
    ('center',          [0.30, -0.30, 13.0, 13.0, 13.0, 13.0, 21.0, 250.0]),
    ('Tb_high_edge',    [0.30, -0.30, 13.0, 13.0, 13.0, 13.0, 21.0, 250.965]),
    ('Tb_low_edge',     [0.30, -0.30, 13.0, 13.0, 13.0, 13.0, 21.0, 249.0]),
    ('Petricca_loweta', [0.30, -0.30, 10.0, 10.5, 10.5, 10.5, 21.0, 250.5]),
    ('cold_stiff',      [0.30, -0.30, 16.0, 16.0, 16.0, 16.0, 23.0, 249.5]),
    ('alpha_low',       [0.10, -0.30, 13.0, 13.0, 13.0, 13.0, 21.0, 250.5]),
    ('alpha_high',      [0.40, -0.30, 13.0, 13.0, 13.0, 13.0, 21.0, 250.5]),
    ('zeta_high',       [0.30,  0.10, 13.0, 13.0, 13.0, 13.0, 21.0, 250.5]),
    ('zeta_low',        [0.30, -1.00, 13.0, 13.0, 13.0, 13.0, 21.0, 250.5]),
]


def main():
    print('Loading structure grid...')
    t0 = time.time()
    grid = build_or_load_structure_grid()
    t_load = time.time() - t0
    print(f'  loaded in {t_load:.1f} s')
    print(f'  Tb_K_grid shape: {np.asarray(grid["Tb_K_grid"]).shape}')
    print(f'  param names:    {PARAM_NAMES}')
    print()

    print(f'{"case":<22s} {"Re(k2)":>10s} {"Im(k2)":>11s} {"|k2|":>8s} '
          f'{"phaseLag_deg":>13s} {"t_ms":>7s}  status')
    print('-' * 92)

    n_pass = n_reject = n_fail = 0
    for name, theta in CASES:
        t0 = time.time()
        try:
            Re_k2, Im_k2, perPhase = forward_model(np.array(theta), grid)
            dt_ms = 1e3 * (time.time() - t0)
            if not np.isfinite(Re_k2):
                status = 'REJECTED (no-ocean safeguard)'
                n_reject += 1
                print(f'{name:<22s} {"nan":>10s} {"nan":>11s} {"nan":>8s} '
                      f'{"nan":>13s} {dt_ms:7.1f}  {status}')
            else:
                k2_mag = np.hypot(Re_k2, Im_k2)
                phase_deg = np.degrees(np.arctan2(-Im_k2, Re_k2))
                status = 'OK'
                n_pass += 1
                print(f'{name:<22s} {Re_k2:10.4f} {Im_k2:11.4e} {k2_mag:8.4f} '
                      f'{phase_deg:13.4f} {dt_ms:7.1f}  {status}')
        except Exception as e:
            dt_ms = 1e3 * (time.time() - t0)
            n_fail += 1
            print(f'{name:<22s} {"ERR":>10s} {"":>11s} {"":>8s} {"":>13s} '
                  f'{dt_ms:7.1f}  FAIL: {type(e).__name__}: {str(e)[:60]}')

    print('-' * 92)
    print(f'pass={n_pass}  rejected={n_reject}  failed={n_fail}  '
          f'(of {len(CASES)} cases)')
    return 0 if n_fail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
