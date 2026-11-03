"""Full PP default runs with Ocean.phaseType lookup vs calc (direct SeaFreeze whichphase). Usage: HOME=~/scratch/pphome python plans/scripts/seafreeze_phase_mode_runs.py Europa Enceladus Ganymede"""
import sys, copy, time, logging, json
sys.path.insert(0, '.')
from importlib import import_module
from PlanetProfile.Main import PlanetProfile
from PlanetProfile.GetConfig import Params as P0
out = {}
for body in sys.argv[1:]:
    for mode in ('lookup', 'calc'):
        logging.disable(logging.WARNING)
        pl = copy.deepcopy(import_module(f'PlanetProfile.Default.{body}.PP{body}').Planet)
        pr = copy.deepcopy(P0)
        pr.CALC_NEW = True; pr.CALC_NEW_GRAVITY = True; pr.CALC_NEW_INDUCT = False
        pr.NO_SAVEFILE = True; pr.SKIP_PLOTS = True; pr.FORCE_EOS_RECALC = True
        pl.Ocean.phaseType = mode
        t = time.perf_counter()
        try:
            pl, pr = PlanetProfile(pl, pr)
            dt = time.perf_counter() - t
            r = dict(t_s=dt, zb_km=float(pl.zb_km), D_ocean=float(pl.D_km), CMR2=float(pl.CMR2mean),
                     k2=[pl.Gravity.k.real, pl.Gravity.k.imag], valid=bool(pl.Do.VALID))
        except Exception as e:
            r = {'error': repr(e)[:160], 't_s': time.perf_counter() - t}
        logging.disable(logging.NOTSET)
        out[f'{body}/{mode}'] = r
        print(body, mode, r, flush=True)
json.dump(out, open('validation_reports/seafreeze_direct_vs_lookup/phase_mode_runs.json', 'w'), indent=1)
