"""Generate the amended joint-Ae anchor configs for the Europa Clipper v2
limits gate (plan Phase 3.4, scientific review 2026-07-13).

Design (why joint, why synthetic):
  A change in Ae(f) moves ALL components {x,y,z} of that frequency together
  (they share the single complex scalar Ae(f); only the fixed excitation
  constants Be_comp differ). Single-channel anchors would sit off the physical
  manifold. So each anchor scales the FIDUCIAL complex Ae(f) to a target
  |Ae| PRESERVING PHASE (isolating the amplitude effect, which is the swept
  quantity), holds Be fixed, and recomputes every kept Bind component of that
  frequency jointly:  Bind_comp = (|Ae|_target * e^{i*phase_fid(f)}) * Be_comp.
  Non-swept frequencies stay at their fiducial Bind values.

  |Ae_synodic| physically spans only ~0.754-0.939 over the induction-supported
  grid, so targets 0.70 and 0.95 are deliberately OFF-manifold. That is valid
  here: the limits gate is a flow-vs-MCMC *consistency* check at each x_obs
  (reject_outside_prior=False), not a claim that each x is realizable by a
  single Tb. Both estimators see the identical (possibly off-manifold) x.

Sweeps (plan-pre-registered):
  (a) synodic  |Ae| in {0.70,0.75,0.80,0.85,0.90,0.95}   (6 anchors)
  (b) synodic 2nd  |Ae| in {0.75,0.83,0.91}  (3, coarse)
  (c) orbital  |Ae| in {0.45,0.62,0.80}      (3, coarse; brackets fiducial 0.62)
  (d) one joint-corner: synodic 0.75 + synodic2nd 0.83 + orbital 0.80

Writes configs to PlanetProfile/Inference/configs/europa_clipper_v2_anchors/.
Each config = v2 config with only the swept Bind_ observables replaced; sigma
(1.5) and all priors/derived/settings unchanged. Emits an anchor manifest
JSON mapping target -> config path + the swept |Ae| for the gate driver.
"""
import os, json, copy
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

V2_CFG = 'PlanetProfile/Inference/configs/europa_seawater_andrade_clipper_v2.json'
OUT_DIR = 'PlanetProfile/Inference/configs/europa_clipper_v2_anchors'
FIDUCIAL_TB = 264.5201
# Kept components per frequency (from the pre-registered pruning; z dropped for
# synodic-2nd and orbital).
KEPT = {
    'synodic':     ['x', 'y', 'z'],
    'synodic 2nd': ['x', 'y'],
    'orbital':     ['x', 'y'],
}
SWEEPS = {
    'synodic':     [0.70, 0.75, 0.80, 0.85, 0.90, 0.95],
    'synodic 2nd': [0.75, 0.83, 0.91],
    'orbital':     [0.45, 0.62, 0.80],
}
CORNER = {'synodic': 0.75, 'synodic 2nd': 0.83, 'orbital': 0.80}

os.makedirs(OUT_DIR, exist_ok=True)

cfg = InferenceConfig.from_json(V2_CFG)
runner = MCMCRunner(cfg)
Tb = np.asarray(runner.structure_data['Tb_K_grid'], float)
i_fid = int(np.argmin(np.abs(Tb - FIDUCIAL_TB)))
Ae_fid = {lbl: complex(runner._ae_grid_cache[i_fid][lbl]) for lbl in KEPT}
be = runner._be_excitation
base = json.load(open(V2_CFG))

def bind_channels(ae_by_label):
    """Bind_<lbl>_<comp>_<part> = (Ae(lbl) * Be_comp).real/.imag for kept comps."""
    out = {}
    for lbl, comps in KEPT.items():
        A = ae_by_label[lbl]
        for c in comps:
            Bind = A * be[lbl][c]
            out[f'Bind_{lbl}_{c}_real'] = round(float(Bind.real), 4)
            out[f'Bind_{lbl}_{c}_imag'] = round(float(Bind.imag), 4)
    return out

def scaled_ae(lbl, target_amp):
    """Scale fiducial complex Ae(lbl) to |Ae|=target, preserving phase."""
    A = Ae_fid[lbl]
    return (target_amp / abs(A)) * A

def write_anchor(name, ae_by_label, meta_extra):
    c = copy.deepcopy(base)
    # Replace ONLY the Bind_ observables; keep sigma from the v2 config.
    new_bind = bind_channels(ae_by_label)
    for k, v in new_bind.items():
        sigma = base['observables'][k][1]
        c['observables'][k] = [v, sigma]
    c['metadata'] = dict(c.get('metadata', {}))
    c['metadata']['anchor_provenance_2026_07_14'] = meta_extra
    path = os.path.join(OUT_DIR, f'{name}.json')
    json.dump(c, open(path, 'w'), indent=2)
    return path

manifest = {'fiducial_tb_K': float(Tb[i_fid]), 'sweeps': {}, 'corner': None,
            'design': 'phase-preserving amplitude scaling of fiducial complex Ae; '
                      'Be fixed; joint over kept comps; non-swept freqs at fiducial'}

for lbl, targets in SWEEPS.items():
    rows = []
    for t in targets:
        ae = dict(Ae_fid)                    # start all frequencies at fiducial
        ae[lbl] = scaled_ae(lbl, t)          # displace only this frequency
        tag = lbl.replace(' ', '') + f'_{t:.2f}'
        meta = (f'Joint-Ae anchor: |Ae_{lbl}|={t:.2f} (phase-preserving scale of '
                f'fiducial Ae_{lbl}={Ae_fid[lbl]:.4f}); other freqs at fiducial. '
                f'reject_outside_prior=False in the gate; off-manifold targets are '
                f'intentional consistency probes.')
        path = write_anchor(f'anchor_{tag}', ae, meta)
        rows.append({'target_amp': t, 'config': path})
    manifest['sweeps'][lbl] = rows

# Joint-corner: displace all three simultaneously.
ae = {lbl: scaled_ae(lbl, CORNER[lbl]) for lbl in KEPT}
meta = (f'Joint-corner anchor: '
        + ', '.join(f'|Ae_{l}|={CORNER[l]:.2f}' for l in KEPT)
        + ' simultaneously (phase-preserving). Stress test of the full manifold corner.')
path = write_anchor('anchor_corner', ae, meta)
manifest['corner'] = {'targets': CORNER, 'config': path}

json.dump(manifest, open(os.path.join(OUT_DIR, 'anchor_manifest.json'), 'w'), indent=2)
print(f"Wrote {sum(len(v) for v in manifest['sweeps'].values()) + 1} anchor configs to {OUT_DIR}")
for lbl, rows in manifest['sweeps'].items():
    print(f"  {lbl}: {[r['target_amp'] for r in rows]}")
print(f"  corner: {CORNER}")

# Sanity: show synodic-x Bind across the sweep (should scale ~linearly in |Ae|).
print("\nsynodic sweep, Bind_synodic_x (real,imag):")
for t in SWEEPS['synodic']:
    ae = dict(Ae_fid); ae['synodic'] = scaled_ae('synodic', t)
    B = ae['synodic'] * be['synodic']['x']
    print(f"  |Ae|={t:.2f}: {B.real:+.2f}{B.imag:+.2f}j")
