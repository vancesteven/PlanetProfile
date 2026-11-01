"""ALMA vs TidalPy gravity-backend baseline on PP default bodies (reviewer
item 6, 2026-11-01). Usage: HOME=~/scratch/pphome python
plans/scripts/gravity_backend_baseline.py Europa Enceladus Titan ..."""
import sys, copy, logging, json
sys.path.insert(0, '.')
from importlib import import_module
from PlanetProfile.Main import PlanetProfile
from PlanetProfile.GetConfig import Params as P0
rows = []
for body in sys.argv[1:]:
    row = {'body': body}
    for backend in ('pyalma', 'tidalpy'):
        logging.disable(logging.WARNING)
        try:
            Planet = copy.deepcopy(import_module(f'PlanetProfile.Default.{body}.PP{body}').Planet)
            Params = copy.deepcopy(P0)
            Params.CALC_NEW = True; Params.CALC_NEW_GRAVITY = True; Params.CALC_NEW_INDUCT = False
            Params.NO_SAVEFILE = True; Params.SKIP_PLOTS = True; Params.QUIET_ALMA = True
            Params.Gravity.backend = backend
            Planet, Params = PlanetProfile(Planet, Params)
            G = Planet.Gravity
            heat = getattr(G, 'tidalpy_Htidal_perPhase_W', None)
            row[backend] = {'backend_used': getattr(G, 'backend_used', None), 'k2': [G.k.real, G.k.imag], 'h2': [G.h.real, G.h.imag],
                            'l2': [G.l.real, G.l.imag], 'libration_m': float(G.libration_m),
                            'heat_W': ({k: float(v) for k, v in heat.items()} if heat else None)}
        except Exception as e:
            row[backend] = {'error': repr(e)[:200]}
        logging.disable(logging.NOTSET)
    rows.append(row)
    a, t = row['pyalma'], row['tidalpy']
    if 'error' in a or 'error' in t:
        print(body, 'ERROR', a.get('error'), t.get('error'), flush=True); continue
    print(f"{body:9s} k2 {a['k2'][0]:.4f}{a['k2'][1]:+.5f}i | {t['k2'][0]:.4f}{t['k2'][1]:+.5f}i   "
          f"h2 {a['h2'][0]:.4f} | {t['h2'][0]:.4f}   lib {a['libration_m']:.1f} | {t['libration_m']:.1f} m  used {t['backend_used']}", flush=True)
json.dump(rows, open(sys.argv[0].replace('.py', '.json') if False else 'validation_reports/gravity_backend_baseline/baseline.json', 'w'), indent=1)
