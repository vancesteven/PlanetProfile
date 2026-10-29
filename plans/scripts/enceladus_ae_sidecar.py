#!/usr/bin/env python
"""Precompute the induction response Ae at every node of the Enceladus (zb, w)
production cache, for DISPLAY ONLY (config metadata.induction_provenance:
no conditioning, no support cut).

Periods are the Moments.py Enceladus set: 'true anomaly' (32.93 h, the
orbital/eccentricity channel, primary per the freeze ruling) and 'synodic'
(15.5592 h, Saturn PPO; non-stationary driver, response only). The cache's
own Texc_hr carries only 'true anomaly', so the periods are passed
explicitly rather than read from each node.

Two core-conductivity curves (freeze ruling, plans/active/enceladus-config-
freeze.md): 'ae_grid' uses PlanetProfile's self-consistent porous-rock sigma
as built; 'ae_grid_core25' imposes 25 S/m on every layer below the ocean
floor (r <= R - zb - D_ocean), the Saur et al. 2024 conductive-core scenario.
The 25 S/m curve is an imposed overlay, NOT a sampled nuisance.

Output uses the MCMCRunner._precompute_ae_grid sidecar schema
({'key': {'labels', 'n', 'cache_size'}, 'ae_grid': {i: {label: complex} |
None}}) so a runner that requests these labels reloads it unchanged.

  python plans/scripts/enceladus_ae_sidecar.py --cache CACHE --out OUT [--workers 10]

Resumable: finished nodes are checkpointed to OUT + '.partial' every 100
nodes; a rerun with the same arguments skips them.
"""
import argparse
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from PlanetProfile.MagneticInduction.Moments import Excitations  # noqa: E402

_STRUCTS = None
_FREQ = None


def _init(cache, freq):
    global _STRUCTS, _FREQ
    with open(cache, 'rb') as f:
        _STRUCTS = pickle.load(f)['structures']
    _FREQ = freq


CORE25_SM = 25.0


def core25(s):
    """Copy of node ``s`` with sigma = 25 S/m in every layer below the ocean
    floor; None if the node has no ocean."""
    import numpy as np
    if not s.get('has_ocean', True) or not s.get('D_ocean_km'):
        return None
    r = np.asarray(s['rSigChange_m'], float)
    sig = np.array(s['sigmaLayers_Sm'], float)
    r_floor = s['R_body_m'] - 1e3 * (s['zb_km_actual'] + s['D_ocean_km'])
    core = r <= r_floor + 1.0
    sig[core] = CORE25_SM
    return {**s, 'sigmaLayers_Sm': sig}


def _one(i):
    from PlanetProfile.Inference.forward_models import forward_model_induction
    s = _STRUCTS[i]
    if not s:
        return i, None, None
    pp = forward_model_induction(s, dict(_FREQ), nn=1, do_parallel=False)
    s25 = core25(s)
    c25 = (forward_model_induction(s25, dict(_FREQ), nn=1, do_parallel=False)
           if s25 is not None else None)
    return i, pp, c25


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=10)
    a = ap.parse_args()
    texc = Excitations.Texc_hr['Enceladus']
    freq = {lab: float(texc[lab]) for lab in ('true anomaly', 'synodic')}
    with open(a.cache, 'rb') as f:
        n = len(pickle.load(f)['structures'])
    t0 = time.time()
    part = a.out + '.partial'
    grid, grid25 = {}, {}
    if os.path.exists(part):
        with open(part, 'rb') as f:
            ck = pickle.load(f)
        if ck.get('n') == n and ck.get('periods_hr') == freq:
            grid, grid25 = ck['ae_grid'], ck['ae_grid_core25']
            print(f'resuming: {len(grid)}/{n} nodes already done', flush=True)

    def checkpoint():
        with open(part + '.tmp', 'wb') as f:
            pickle.dump({'n': n, 'periods_hr': freq, 'ae_grid': grid,
                         'ae_grid_core25': grid25}, f)
        os.replace(part + '.tmp', part)

    todo = [i for i in range(n) if i not in grid]
    with ProcessPoolExecutor(a.workers, initializer=_init,
                             initargs=(a.cache, freq)) as ex:
        for k, (i, r, r25) in enumerate(ex.map(_one, todo, chunksize=4)):
            grid[i] = r
            grid25[i] = r25
            if (k + 1) % 100 == 0 or k + 1 == len(todo):
                checkpoint()
                print(f'{len(grid)}/{n} nodes, {time.time() - t0:.0f} s',
                      flush=True)
    payload = {'key': {'labels': tuple(sorted(freq)), 'n': n,
                       'cache_size': os.path.getsize(a.cache)},
               'ae_grid': grid, 'ae_grid_core25': grid25,
               'core25_Sm': CORE25_SM, 'periods_hr': freq,
               'note': 'display only; enceladus_ae_sidecar.py'}
    tmp = a.out + '.tmp'
    with open(tmp, 'wb') as f:
        pickle.dump(payload, f)
    os.replace(tmp, a.out)
    if os.path.exists(part):
        os.remove(part)
    print('wrote', a.out, 'non-None', sum(v is not None for v in grid.values()))


if __name__ == '__main__':
    main()
