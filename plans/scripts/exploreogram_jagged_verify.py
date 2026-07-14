"""Drive a real 5x5 Europa exploreogram through the page (user preset:
sigmaMean/salinity x, D_km/ice-thickness y, induction z, Bx, log axes),
matplotlib path, and report where figures were saved."""
import glob
import os
import sys
import time

sys.path.insert(0, '/Users/svance/ppgenai/PlanetProfileApp')
sys.path.insert(0, '/Users/svance/ppgenai')
os.chdir('/Users/svance/ppgenai/PlanetProfileApp')

from streamlit.testing.v1 import AppTest

t0 = time.time()
at = AppTest.from_file('/Users/svance/ppgenai/PlanetProfileApp/pages/Exploreogram.py',
                       default_timeout=2400)
at.run()
assert not at.exception, [str(e.value) for e in at.exception]

# Small grid for the integration check
for ni in at.number_input:
    if ni.key == 'nx_input':
        ni.set_value(5)
    if ni.key == 'ny_input':
        ni.set_value(5)
    if ni.key == 'x_min':
        ni.set_value(1.0)
    if ni.key == 'x_max':
        ni.set_value(100.0)
at.run()

run_btns = [b for b in at.button if 'Run Exploreogram' in (b.label or '')]
print('buttons:', [b.label for b in at.button][:6])
assert run_btns, 'Run button not found'
run_btns[0].click().run()
print('elapsed: %.1f min' % ((time.time() - t0) / 60))
print('exceptions:', [str(e.value)[:200] for e in at.exception])
print('errors:', [str(e.value)[:200] for e in at.error][:4])
res = at.session_state['explore_results'] if 'explore_results' in at.session_state else None
print('results stored:', res is not None)

# Find figures written in the last 40 minutes
cands = []
for root in ('/Users/svance/ppgenai/figures', '/Users/svance/ppgenai/PlanetProfileApp',
             '/Users/svance/ppgenai/Europa'):
    for pat in ('**/*.png', '**/*.pdf'):
        for f in glob.glob(os.path.join(root, pat), recursive=True):
            if time.time() - os.path.getmtime(f) < 2400:
                cands.append(f)
print('recent figures:', cands[:10])
