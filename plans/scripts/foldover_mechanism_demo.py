"""Before/after render of the fold-over guard on a synthetic grid mimicking
the reported artifact: y = derived quantity with non-monotone wiggle that
grows with x (high salinity), z smooth."""
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/Users/svance/ppgenai')
from PlanetProfile.Plotting.ExplorationPlots import _sort_grid_axes

rng = np.random.default_rng(3)
nx, ny = 24, 24
sal = np.linspace(1, 100, nx)                      # x: salinity-like
ybase = np.linspace(20, 120, ny)                   # y: ocean-thickness-like
x = np.tile(sal[:, None], (1, ny))
# wiggle amplitude grows with salinity -> local non-monotonicity at high x
wiggle = (sal[:, None] / 100.0) ** 2 * 12.0 * np.sin(np.arange(ny) * 2.2)
y = np.tile(ybase[None, :], (nx, 1)) + wiggle + rng.normal(0, 0.5, (nx, ny))
z = np.sqrt(x) * 3 + (y / 20.0) ** 1.5             # smooth field

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
axes[0].pcolormesh(x, y, z, shading='auto', cmap='viridis')
axes[0].set_title('BEFORE: raw non-monotone y (folded quads)')
xs, ys, zs, _ = _sort_grid_axes(x, y, z, None)
axes[1].pcolormesh(xs, ys, zs, shading='auto', cmap='viridis')
axes[1].set_title('AFTER: _sort_grid_axes fold-over guard')
for a in axes:
    a.set_xlabel('salinity-like x')
axes[0].set_ylabel('derived y')
fig.tight_layout()
out = ('/private/tmp/claude-501/-Users-svance-ppgenai/'
       '9777aafb-f62d-4c15-985a-f9e7f3898c20/scratchpad/foldover_demo.png')
fig.savefig(out, dpi=110)
print('WROTE', out)
