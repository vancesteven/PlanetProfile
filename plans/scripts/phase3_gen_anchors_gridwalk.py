"""CORRECTED Europa Clipper v2 limits-gate anchors: single-Tb grid-walk.

Scientific review (2026-07-14, scientific-reviewer) REJECTED the independent
per-frequency phase-preserving-scaling design because in this forward model
Ae is a function of Tb ALONE: the single complex scalar Ae(f, Tb) for each
excitation frequency is set by one Tb via the cached conductivity profile.
Therefore:
  - |Ae| and arg(Ae) are locked together along the physical Ae(Tb) locus
    (synodic Delta-arg ~21 deg, orbital ~45 deg across support); fixed-phase
    radial scaling leaves the manifold everywhere but the fiducial.
  - The three frequencies cannot be moved independently; pinning non-swept
    frequencies at fiducial while moving one is internally inconsistent.
  - The reference MCMC samples Tb and applies a HARD |Ae_synodic|>0.7 cut
    (reject_outside_prior is SBI-only), so off-manifold anchors rail to
    boundary garbage and are useless W1 references.

CORRECT construction (the reviewer's remedy): walk real grid Tb points across
the induction-supported range. At each anchor Tb, set ALL 14 Bind channels
from that ONE Tb's actual complex Ae(f, Tb) * Be_comp(f). Every anchor is on
the physical manifold, realizable by a single Tb, and both estimators operate
on their valid domain. This is the PRIMARY W1 calibration gate.

Anchor Tb points: chosen on the real grid to span |Ae_synodic| across its
achievable support range [~0.754, ~0.939] with good coverage, including the
fiducial neighborhood. Non-swept? -> nothing is "non-swept": all frequencies
co-vary, which is the physically correct joint sweep.
"""
import os, json, copy
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

V2_CFG = 'PlanetProfile/Inference/configs/europa_seawater_andrade_clipper_v2.json'
OUT_DIR = 'PlanetProfile/Inference/configs/europa_clipper_v2_anchors_gridwalk'
KEPT = {'synodic': ['x', 'y', 'z'], 'synodic 2nd': ['x', 'y'], 'orbital': ['x', 'y']}
# Target grid Tb points (K) spanning the induction-supported range. These are
# REQUESTED values; the generator snaps each to the nearest actual grid point.
TB_TARGETS = [261.5, 262.5, 263.5, 264.5, 266.0, 268.0, 270.0, 271.0]

os.makedirs(OUT_DIR, exist_ok=True)
cfg = InferenceConfig.from_json(V2_CFG)
runner = MCMCRunner(cfg)
Tb = np.asarray(runner.structure_data['Tb_K_grid'], float)
be = runner._be_excitation
base = json.load(open(V2_CFG))

def bind_channels_at(idx):
    Ae = runner._ae_grid_cache[idx]
    out = {}
    for lbl, comps in KEPT.items():
        A = complex(Ae[lbl])
        for c in comps:
            B = A * be[lbl][c]
            out[f'Bind_{lbl}_{c}_real'] = round(float(B.real), 4)
            out[f'Bind_{lbl}_{c}_imag'] = round(float(B.imag), 4)
    return out, {lbl: abs(complex(Ae[lbl])) for lbl in KEPT}

manifest = {'design': 'single-Tb grid-walk; all frequencies co-vary from one '
                      'physical Ae(Tb); PRIMARY W1 calibration gate',
            'anchors': []}
used_idx = set()
for tb_req in TB_TARGETS:
    idx = int(np.argmin(np.abs(Tb - tb_req)))
    if idx in used_idx:
        continue
    used_idx.add(idx)
    tb_act = float(Tb[idx])
    bind, amps = bind_channels_at(idx)
    c = copy.deepcopy(base)
    for k, v in bind.items():
        c['observables'][k] = [v, base['observables'][k][1]]
    c['metadata'] = dict(c.get('metadata', {}))
    c['metadata']['anchor_provenance_2026_07_14'] = (
        f'Grid-walk anchor at Tb={tb_act:.3f} K (grid idx {idx}). All 14 Bind '
        f'channels set from this single Tb actual complex Ae(f)*Be_comp; on the '
        f'physical manifold. |Ae|: synodic={amps["synodic"]:.4f}, '
        f'synodic2nd={amps["synodic 2nd"]:.4f}, orbital={amps["orbital"]:.4f}. '
        f'PRIMARY limits/W1 calibration gate (reviewer-corrected design).')
    tag = f'tb{tb_act:.2f}'.replace('.', 'p')
    path = os.path.join(OUT_DIR, f'anchor_{tag}.json')
    json.dump(c, open(path, 'w'), indent=2)
    manifest['anchors'].append({
        'tb_K': tb_act, 'grid_idx': idx, 'config': path,
        'amp_synodic': round(amps['synodic'], 4),
        'amp_synodic2nd': round(amps['synodic 2nd'], 4),
        'amp_orbital': round(amps['orbital'], 4),
    })

json.dump(manifest, open(os.path.join(OUT_DIR, 'anchor_manifest.json'), 'w'), indent=2)
print(f"Wrote {len(manifest['anchors'])} grid-walk anchor configs to {OUT_DIR}")
print("%-10s %-6s %-10s %-12s %-10s" % ('Tb_K', 'idx', '|Ae_syn|', '|Ae_2nd|', '|Ae_orb|'))
for a in manifest['anchors']:
    print("%-10.3f %-6d %-10.4f %-12.4f %-10.4f" % (
        a['tb_K'], a['grid_idx'], a['amp_synodic'], a['amp_synodic2nd'], a['amp_orbital']))
