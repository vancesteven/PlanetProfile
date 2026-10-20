#!/usr/bin/env python
"""READ-ONLY acceptance check of the merged Enceladus (zb, w) production cache.

WHY THIS EXISTS. The production build (plans/scripts/
enceladus_zbw_production_build.py, or the SLURM shards merged by
enceladus_zbw_merge_shards.py) enforces most gates INSIDE the builder and the
merge script reports gate figures at merge time. Neither is an independent
audit of the delivered artifact. This script re-derives every gate from the
saved cache pickle alone (plus the committed config), never re-runs or relaxes
anything, never tunes a threshold, and writes its verdicts next to the cache
as `<cache>.acceptance.json`. Exit status is nonzero iff any gate is FAIL.
UNCHECKABLE (data not present in the artifact) is reported, not counted as
PASS, and does not by itself change the exit status.

Usage:
  HOME=~/scratch/pphome python plans/scripts/enceladus_zbw_acceptance.py \
      --cache ~/scratch/enceladus_zbw_production_cache.pkl \
      [--shard-dir ~/scratch/enc_zbw_shards]

GATE SOURCES (quoted verbatim).

A) plans/MACHINE-B-HANDOFF.md section 0.31 item 4 ("Acceptance (in addition
   to the six section 0.29 gates)"):
   A1 "exactly 0 C2 violations;"
   A2 "(6.0, 2.8943) built with zb_actual 6.0086;"
   A3 "(12.0, 4.9239) None with the unreachable reason (0.067 km below
       its floor, outside the band);"
   A4 "built/None set otherwise identical to the failed r7 run's log."
      "Report any difference rather than reconciling it."
   A4 is UNCHECKABLE: the failed r7 run raised before _save_cache_atomic, so
   no machine-readable built/None record of it was committed (only the prose
   in validation_reports/enceladus_isostasy/r7_build_c2_false_positive.md,
   which names the single offending node). This script dumps THIS cache's
   built/None set to `<cache>.built_none_set.json` for manual comparison
   against the failed run's console log.

B) validation_reports/enceladus_isostasy/r7_ADJUDICATION.md, "MACHINE B BUILD
   INSTRUCTIONS": "ACCEPTANCE: 0 placement rejects, max residual <0.125
   (expect <=0.035); None set == C1 prediction within one cell (any outside
   = halt); re-measured MoI floor at production scale; mass invariant on
   every node; frozen I-F1..I-F4 + I-F6 must NOT fire; f_excluded recorded."
   The adjudication lists these as a run-on sentence; the six numbered gates
   follow enceladus_zbw_merge_shards.py's grouping (its acceptance_gates keys):
   G1 0 placement rejects            (n_zb_placement_rejected == 0)
   G2 max residual < 0.125 km        (expectation <=0.035 reported, not gated)
   G3 None set == C1 prediction within one cell
   G4 re-measured MoI floor at production scale (ocean_moi_window)
   G5 mass invariant on every node; frozen I-F1..I-F4 + I-F6 must NOT fire
   G6 f_excluded recorded (C3 normalization input)
   (The merge script's own gate 5 is frozen-invariants and its gate 3/6 are
   unconditional passes; here each is re-derived from the node data.)

C) plans/MACHINE-B-HANDOFF.md section 0.30 ("INCLUSIVE -- 88 nodes, ending at
   45.0 km"), r7 instructions (40 w nodes, 3520 ocean nodes, frozen axis
   arange(46.5, 65.8001, 0.5) = 39 nodes) and "zb_tol_km=0.125 EXPLICIT".
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

CONFIG_PATH = REPO / ('PlanetProfile/Inference/configs/'
                      'enceladus_cassini_isostasy_7D.json')

# Constants taken from the ratified rulings / builder, not tuned here.
ZB_TOL_KM = 0.125                    # r7 "zb_tol_km=0.125 EXPLICIT"
RESIDUAL_EXPECTATION_KM = 0.035      # r7 "(expect <=0.035)", informational
FROZEN_MASS_TOL = 1e-6               # builder I-F1 default / r7 mass_tol
FROZEN_RHO_CLOSURE_TOL = 12.0        # I-F3, kg/m^3
FROZEN_ZB_TOL_KM = 0.25              # I-F4
OCEAN_MASS_POSITIVE_EPS = 1e-9       # cache_builder._FROZEN_IF2_POSITIVE_EPS:
#   residual is one-signed by design (Silicates.py:118); positive = bug class
N_ZB, N_W, N_OCEAN, N_FROZEN = 88, 40, 3520, 39
ZB_END_KM = 45.0

RESULTS = []


def record(gate, status, value, detail='', source=''):
    assert status in ('PASS', 'FAIL', 'UNCHECKABLE')
    RESULTS.append({'gate': gate, 'status': status, 'value': value,
                    'detail': detail, 'source': source})
    print(f'[{status:11s}] {gate}: {value}' + (f'  -- {detail}' if detail else ''))


def zb_min_interp(w_ppt, table):
    """Same math as cache_builder._reachability_zb_min_interp (log-linear in
    w, clamped), re-stated here so the audit does not import the builder."""
    rows = sorted(table, key=lambda r: float(r['w_ppt']))
    log_w = np.log10(np.asarray([float(r['w_ppt']) for r in rows]))
    zb_min = np.asarray([float(r['zb_min_km']) for r in rows])
    return float(np.interp(np.log10(max(float(w_ppt), 1e-300)), log_w, zb_min))


def find_node(zb, w, zb_t, w_t):
    iz = np.where(np.isclose(zb, zb_t, atol=1e-9))[0]
    iw = np.where(np.abs(w - w_t) < 5e-5)[0]   # table stores w to 4 dp
    if iz.size != 1 or iw.size != 1:
        return None
    return int(iz[0]), int(iw[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--cache', required=True)
    ap.add_argument('--shard-dir', default=None)
    args = ap.parse_args()

    cache_path = Path(args.cache).expanduser()
    with open(cache_path, 'rb') as f:
        cache = pickle.load(f)
    with open(CONFIG_PATH) as f:
        cfg = json.load(f)
    spec = cfg['metadata']['structure_cache_spec']

    print(f'Cache: {cache_path}\nschema_version={cache.get("schema_version")}\n')

    if cache.get('QUARANTINED_C2_FAILED'):
        record('cache_not_quarantined', 'FAIL', True,
               'cache is a QUARANTINED_C2_FAILED dump; never consume it')

    zb = np.asarray(cache['zb_km_grid'], dtype=float)
    w = np.asarray(cache['wOcean_ppt_grid'], dtype=float)
    structs = cache['structures']
    reasons = cache.get('ocean_none_reasons') or [None] * len(structs)
    n_zb, n_w = zb.size, w.size
    built = np.array([s is not None for s in structs])

    # ---------------- Section 0.30 / axis & tolerance checks --------------
    record('0.30_zb_axis_88_nodes_end_45.0', 'PASS' if (
        n_zb == N_ZB and abs(zb[-1] - ZB_END_KM) < 1e-9 and zb[0] == 5.0
        and np.all(np.diff(zb) > 0)) else 'FAIL',
        f'n_zb={n_zb}, first={zb[0]}, last={zb[-1]}',
        'inclusive axis ending at 45.0 km', 'MACHINE-B-HANDOFF 0.30')
    record('r7_w_axis_40_nodes', 'PASS' if (
        n_w == N_W and np.array_equal(
            w, 10 ** np.linspace(-1.0, 2.0, 40))) else 'FAIL',
        f'n_w={n_w}', 'must be bit-identical to 10**linspace(-1,2,40)',
        'r7_ADJUDICATION build instructions')
    record('0.30_ocean_nodes_3520', 'PASS' if (
        len(structs) == N_OCEAN == n_zb * n_w) else 'FAIL',
        f'len(structures)={len(structs)}', '', 'MACHINE-B-HANDOFF 0.30')
    fz = cache.get('frozen_zb_km_grid')
    frozen = cache.get('frozen_structures')
    record('r7_frozen_nodes_39_axis', 'PASS' if (
        fz is not None and len(fz) == N_FROZEN and frozen is not None
        and len(frozen) == N_FROZEN
        and np.array_equal(np.asarray(fz), np.arange(46.5, 65.8001, 0.5)))
        else 'FAIL', f'n_frozen={None if fz is None else len(fz)}',
        'arange(46.5, 65.8001, 0.5)', 'r7_ADJUDICATION build instructions')
    ztol = cache.get('zb_tol_km')
    record('r7_zb_tol_km_0.125', 'PASS' if (
        ztol is not None and abs(float(ztol) - ZB_TOL_KM) < 1e-12)
        else 'FAIL', f'zb_tol_km={ztol}', '', 'r7 "zb_tol_km=0.125 EXPLICIT"')
    cfg_n = (spec['zb_km_grid']['n_total'], spec['log10_w_grid']['n'])
    record('config_declares_88x40', 'PASS' if cfg_n == (N_ZB, N_W) else 'FAIL',
           f'config n_total/n_w={cfg_n}', '', 'config structure_cache_spec')

    # ---------------- A) section 0.31 acceptance points -------------------
    violations = []
    table = cache.get('ocean_reachability_restriction_table')
    band = cache.get('reachability_band_km')
    pred = cache.get('reachability_predicate')
    if not table or band is None:
        record('A1_zero_C2_violations', 'UNCHECKABLE', None,
               'cache lacks ocean_reachability_restriction_table / '
               'reachability_band_km (C2 did not run, or pre-0.31 cache)',
               'cache_builder.py build_zbw_grid_cache return dict')
    else:
        cfg_table = spec['ocean_reachability_restriction']['table']
        record('A0_cache_table_equals_config_table',
               'PASS' if table == cfg_table and len(table) == N_W else 'FAIL',
               f'rows={len(table)}', 'cache table must equal committed config',
               'config structure_cache_spec.ocean_reachability_restriction')
        record('A0_band_is_0.4_zb_tol', 'PASS' if (
            abs(float(band) - 0.4 * ZB_TOL_KM) < 1e-12) else 'FAIL',
            f'reachability_band_km={band}; predicate={pred!r}',
            'band = 0.4*0.125 = 0.05 km', 'MACHINE-B-HANDOFF 0.31 item 2')
        violations, excluded = [], np.zeros(len(structs), dtype=bool)
        for iz in range(n_zb):
            for iw in range(n_w):
                k = iz * n_w + iw
                zmin = zb_min_interp(w[iw], table)
                ex = zb[iz] < zmin - float(band)
                excluded[k] = ex
                if ex and built[k]:
                    violations.append(f'BUILT but excluded (zb={zb[iz]}, '
                                      f'w={w[iw]:.4f}, floor={zmin:.4f})')
                elif (not ex and not built[k]
                      and reasons[k] != 'no_ocean_frozen_edge'):
                    violations.append(f'None({reasons[k]}) but reachable '
                                      f'(zb={zb[iz]}, w={w[iw]:.4f}, '
                                      f'floor={zmin:.4f})')
        record('A1_zero_C2_violations', 'PASS' if not violations else 'FAIL',
               len(violations), '; '.join(violations[:5]),
               'MACHINE-B-HANDOFF 0.31 item 4 bullet 1')

        # A2
        n62 = find_node(zb, w, 6.0, 2.8943)
        if n62 is None:
            record('A2_node_6.0_2.8943_built_zb_actual_6.0086', 'FAIL', None,
                   'node not found on axes')
        else:
            s = structs[n62[0] * n_w + n62[1]]
            za = None if s is None else float(s.get('zb_km_actual', np.nan))
            record('A2_node_6.0_2.8943_built_zb_actual_6.0086',
                   'PASS' if (s is not None and za is not None
                              and round(za, 4) == 6.0086) else 'FAIL',
                   f'built={s is not None}, zb_actual={za}',
                   'compared at 4 dp (table precision)',
                   'MACHINE-B-HANDOFF 0.31 item 4 bullet 2')
        # A3
        n124 = find_node(zb, w, 12.0, 4.9239)
        if n124 is None:
            record('A3_node_12.0_4.9239_None_unreachable', 'FAIL', None,
                   'node not found on axes')
        else:
            k = n124[0] * n_w + n124[1]
            zmin = zb_min_interp(w[n124[1]], table)
            deficit = zmin - 12.0
            ok = (structs[k] is None and bool(excluded[k])
                  and reasons[k] == 'build_exception'
                  and deficit > float(band))
            record('A3_node_12.0_4.9239_None_unreachable',
                   'PASS' if ok else 'FAIL',
                   f'None={structs[k] is None}, reason={reasons[k]!r}, '
                   f'floor={zmin:.4f}, below-floor={deficit:.4f} km '
                   f'(handoff: ~0.067; band {band})',
                   "reason 'build_exception' is the only reason the builder "
                   "assigns to a brentq same-sign-bracket failure",
                   'MACHINE-B-HANDOFF 0.31 item 4 bullet 3')

    # A4: no machine-readable failed-run record exists.
    bn = {'n_zb': n_zb, 'n_w': n_w,
          'zb_km_grid': zb.tolist(), 'wOcean_ppt_grid': w.tolist(),
          'none_flat_indices': [int(i) for i in np.where(~built)[0]],
          'none_reasons': {str(int(i)): reasons[i] for i in np.where(~built)[0]},
          'flat_index': 'i_zb * n_w + i_w'}
    bn_path = Path(str(cache_path) + '.built_none_set.json')
    with open(bn_path, 'w') as f:
        json.dump(bn, f)
    record('A4_built_none_set_identical_to_failed_r7_log', 'UNCHECKABLE',
           f'{int((~built).sum())} None / {int(built.sum())} built',
           'no machine-readable record of the failed r7 run exists (it died '
           'before saving; only prose in r7_build_c2_false_positive.md and '
           'the console log on the build machine). This cache set was written '
           f'to {bn_path} -- diff it against that log by hand and report any '
           'difference without reconciling.',
           'MACHINE-B-HANDOFF 0.31 item 4 bullet 4')

    # ---------------- B) six section 0.29 gates ---------------------------
    nrej = int(cache.get('n_zb_placement_rejected', -1))
    resid = [abs(float(s['zb_residual_km'])) for s in structs if s is not None]
    record('G1_zero_placement_rejects', 'PASS' if nrej == 0 else 'FAIL', nrej,
           "'0 placement rejects'", 'r7_ADJUDICATION ACCEPTANCE')
    nrej_reason = sum(1 for r in reasons if r == 'zb_placement_reject')
    record('G1b_placement_reject_reasons_zero',
           'PASS' if nrej_reason == 0 else 'FAIL', nrej_reason,
           'count of ocean_none_reasons == zb_placement_reject')
    mx = max(resid) if resid else float('nan')
    record('G2_max_residual_lt_0.125', 'PASS' if (resid and mx < ZB_TOL_KM)
           else 'FAIL', f'max|residual|={mx:.5f} km',
           f"'max residual <0.125 (expect <=0.035)'; expectation "
           f"{'met' if mx <= RESIDUAL_EXPECTATION_KM else 'NOT met (informational)'}",
           'r7_ADJUDICATION ACCEPTANCE')

    edge = [(zb[k // n_w]) for k in range(len(structs))
            if not built[k] and reasons[k] == 'no_ocean_frozen_edge']
    if not table:
        record('G3_none_set_equals_C1_prediction', 'UNCHECKABLE', None,
               'no table in cache')
    else:
        # Exact predicate mismatch count is A1's. "Within one cell": report
        # any mismatch that would be absorbed by shifting to a neighbouring zb
        # node, but the exact result stays the gate (never tuned).
        n_exact = len(violations)
        record('G3_none_set_equals_C1_prediction',
               'PASS' if n_exact == 0 else 'FAIL',
               f'{n_exact} nodes outside prediction (exact predicate); '
               f'{int((~built).sum())} total None, of which '
               f'{len(edge)} are no_ocean_frozen_edge '
               f'(zb {min(edge) if edge else None}-{max(edge) if edge else None} km)',
               "'None set == C1 prediction within one cell (any outside = "
               "halt)'; edge Nones are the separately declared ~45 km edge",
               'r7_ADJUDICATION ACCEPTANCE')

    win = cache.get('ocean_moi_window')
    if not win:
        record('G4_moi_floor_remeasured', 'UNCHECKABLE', None,
               'ocean_moi_window absent')
    else:
        cm = win.get('Cmeasured')
        devs = [abs(float(s['CMR2']) - cm) for s in structs
                if s is not None and s.get('CMR2') is not None
                and np.isfinite(float(s['CMR2']))] if cm is not None else []
        rec_max = win.get('max_abs_cmr2_deviation')
        re_max = max(devs) if devs else float('nan')
        half_u, half_l = win.get('CuncertaintyUpper'), win.get('CuncertaintyLower')
        ok = (rec_max is not None and np.isfinite(rec_max)
              and abs(rec_max - re_max) < 1e-12
              and half_u == 0.08 and half_l == 0.08
              and win.get('margin', -1) >= 0
              and win.get('n_nodes_measured') == int(built.sum())
              and str(win.get('source', '')).startswith('config metadata'))
        record('G4_moi_floor_remeasured', 'PASS' if ok else 'FAIL',
               f"recorded max|dC|={rec_max}, recomputed={re_max}, "
               f"half-width={half_u}/{half_l}, margin={win.get('margin')}, "
               f"n_measured={win.get('n_nodes_measured')} "
               f"(built {int(built.sum())}), source={win.get('source')!r}",
               'recorded == recomputed from node CMR2; Cuncertainty 0.08; '
               'margin >= 0; provenance = config',
               'r7_ADJUDICATION C-provenance chain / ACCEPTANCE')

    # G5 mass invariant on every node
    om = [float(s.get('mass_residual_frac', np.nan)) for s in structs
          if s is not None]
    o_ok = bool(om) and np.all(np.isfinite(om)) and max(om) <= OCEAN_MASS_POSITIVE_EPS
    record('G5a_ocean_mass_invariant_every_node', 'PASS' if o_ok else 'FAIL',
           f'n={len(om)}, min={min(om) if om else None}, '
           f'max={max(om) if om else None}',
           f'finite on every built node and not positive beyond '
           f'{OCEAN_MASS_POSITIVE_EPS:g} (solver residual is one-signed by '
           'design; no ratified numeric lower bound, so min is reported only)',
           'cache_builder.py build_single_structure docstring / r7 ACCEPTANCE')

    fbs = cache.get('frozen_build_spec') or {}
    fnodes = [s for s in (frozen or []) if s is not None]
    fm = [abs(float(s.get('mass_residual_frac', np.nan))) for s in fnodes]
    record('G5b_frozen_mass_invariant_I-F1_every_node', 'PASS' if (
        fnodes and len(fnodes) == N_FROZEN and np.all(np.isfinite(fm))
        and max(fm) <= FROZEN_MASS_TOL) else 'FAIL',
        f'n={len(fnodes)}, max|res|={max(fm) if fm else None}',
        f'<= {FROZEN_MASS_TOL:g}', 'cache_builder.py I-F1')
    rej = fbs.get('rejects_by_invariant')
    record('G6a_frozen_no_rejects_I-F1..I-F4', 'PASS' if (
        rej is not None and sum(rej.values()) == 0
        and fbs.get('n_rejected') == 0 and fbs.get('n_built') == N_FROZEN
        and cache.get('n_frozen_if2_rejected') == 0) else 'FAIL',
        f"rejects_by_invariant={rej}, n_built={fbs.get('n_built')}, "
        f"n_rejected={fbs.get('n_rejected')}, "
        f"n_frozen_if2_rejected={cache.get('n_frozen_if2_rejected')}",
        "'frozen I-F1..I-F4 ... must NOT fire'", 'r7_ADJUDICATION ACCEPTANCE')
    clo = [abs(float(s.get('rho_rock_closure_residual_kgm3', np.nan)))
           for s in fnodes]
    fzr = [abs(float(s.get('zb_residual_km', np.nan))) for s in fnodes]
    nliq = [int(s.get('n_liquid_layers', -1)) for s in fnodes]
    ok = (fnodes and np.all(np.isfinite(clo)) and max(clo) <= FROZEN_RHO_CLOSURE_TOL
          and np.all(np.isfinite(fzr)) and max(fzr) < FROZEN_ZB_TOL_KM
          and all(n == 0 for n in nliq))
    record('G6b_frozen_nodes_recheck_I-F2liquid_I-F3_I-F4', 'PASS' if ok else 'FAIL',
           f'max closure resid={max(clo) if clo else None} kg/m3 '
           f'(<= {FROZEN_RHO_CLOSURE_TOL}), max|zb resid|='
           f'{max(fzr) if fzr else None} km (< {FROZEN_ZB_TOL_KM}), '
           f'liquid layers all 0={all(n == 0 for n in nliq)}',
           'recomputed from stored per-node fields',
           'cache_builder.py I-F2_liquid/I-F3/I-F4')
    # I-F6: not every node inside the template MoI window
    win6 = fbs.get('I-F6_window')
    cm6 = [float(s.get('cmr2_derived', np.nan)) for s in fnodes]
    if win6 and cm6:
        inside = np.abs(np.asarray(cm6) - win6[0]) <= win6[1]
        record('G6c_frozen_I-F6_non_conditioning', 'PASS' if not np.all(inside)
               else 'FAIL', f'{int(inside.sum())}/{len(cm6)} nodes inside '
               f'{win6[0]} +/- {win6[1]}; span {min(cm6):.5f}-{max(cm6):.5f}',
               'build fails iff ALL nodes inside window',
               'cache_builder.py I-F6')
    else:
        record('G6c_frozen_I-F6_non_conditioning', 'UNCHECKABLE', None,
               'frozen_build_spec.I-F6_window or cmr2_derived missing')

    # G6 (merge-script gate 6) f_excluded recorded
    n_none = int((~built).sum())
    f_ex = n_none / len(structs)
    evp = cfg['metadata'].get('branch_model', {}).get('evidence_protocol', {})
    has_c3 = any('reachability_normalization' in k for k in evp) or any(
        'reachability_normalization' in k
        for k in spec.get('ocean_reachability_restriction', {}))
    rep = REPO / 'validation_reports/enceladus_isostasy/r7_production_build_report.json'
    rep_f = None
    if rep.exists():
        try:
            rep_f = json.load(open(rep))['acceptance_gates'][
                '6_f_excluded_recorded']['f_excluded']
        except Exception:
            rep_f = None
    if rep_f is not None and abs(rep_f - f_ex) > 1e-12:
        record('G7_f_excluded_recorded', 'FAIL', f'cache f_excluded={f_ex:.6f}',
               f'merge report records {rep_f}', 'r7_ADJUDICATION C3')
    else:
        record('G7_f_excluded_recorded', 'PASS' if has_c3 else 'FAIL',
               f'f_excluded = {n_none}/{len(structs)} = {f_ex:.6f} '
               f'(1-f = {1 - f_ex:.6f}; C3 shift +log(1/(1-f)) = '
               f'{-np.log1p(-f_ex):.6f}); merge report '
               f'{"agrees" if rep_f is not None else "not present"}',
               'computed from the cache None mask (the cache stores no '
               'f_excluded field); C3 config validity-condition text present: '
               f'{has_c3}. NOTE: n_none includes any no_ocean_frozen_edge '
               'Nones -- C3 wants reachability and MoI/edge exclusions '
               'reported separately; see G3 for the edge count.',
               'r7_ADJUDICATION C3 / ACCEPTANCE')

    # ---------------- optional shard cross-check ---------------------------
    if args.shard_dir:
        sd = Path(args.shard_dir).expanduser()
        metas = sorted(sd.glob('shard_*_of_*.meta.json'))
        quar = list(sd.glob('*preC2_QUARANTINED*')) + list(
            cache_path.parent.glob('*preC2_QUARANTINED*'))
        try:
            all_st, all_n, ok = [], [], bool(metas)
            for m in metas:
                meta = json.load(open(m))
                p = Path(meta['cache_path'])
                if not p.exists():
                    p = m.parent / m.name.replace('.meta.json', '.pkl')
                sc = pickle.load(open(p, 'rb'))
                all_st.extend(s is not None for s in sc['structures'])
            same = (np.array(all_st) == built).all() if len(all_st) == len(
                built) else False
            record('S_shards_complete_and_match_merged',
                   'PASS' if (metas and same and not quar) else 'FAIL',
                   f'{len(metas)} shards, built/None mask equal={bool(same)}, '
                   f'quarantine files={[str(q) for q in quar]}',
                   'shard concatenation vs merged structures')
        except Exception as exc:  # report, never crash the audit
            record('S_shards_complete_and_match_merged', 'FAIL', None,
                   f'{type(exc).__name__}: {exc}')
    else:
        record('S_shards_complete_and_match_merged', 'UNCHECKABLE', None,
               '--shard-dir not supplied')

    # ---------------- summary ----------------------------------------------
    counts = {s: sum(r['status'] == s for r in RESULTS)
              for s in ('PASS', 'FAIL', 'UNCHECKABLE')}
    out = Path(str(cache_path) + '.acceptance.json')
    with open(out, 'w') as f:
        json.dump({'cache': str(cache_path), 'counts': counts,
                   'results': RESULTS}, f, indent=2, default=str)
    print(f'\nSUMMARY {counts}  -> {out}')
    print('Gates are interpreted, never tuned. Report FAIL/UNCHECKABLE to '
          'Machine A as-is.')
    sys.exit(1 if counts['FAIL'] else 0)


if __name__ == '__main__':
    main()
