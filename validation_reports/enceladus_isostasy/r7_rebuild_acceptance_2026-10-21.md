# Enceladus r7 production rebuild — merge + acceptance (2026-10-21)

- **Build:** third local launch on the Mac, 11 shards. BUILD_COMMIT 634ea0cd, which comes after the §0.31 C2 fix cf0adfb7. DONE 2026-10-20 20:48 BST. No Traceback, no C2 quarantine dump.
- **Merge and acceptance:** `plans/cluster/enceladus_zbw_finish.sh` at c81e7398, run 2026-10-21 11:22 BST.
- **Merged cache:** `~/scratch/enceladus_zbw_production_cache.pkl` on the Mac. It is not yet in the repo and waits on adjudication.
- **Mac log:** `~/scratch/enceladus_zbw_finish.20261021_112221.log`.

## Summary as run: 23 PASS, 0 FAIL, 1 UNCHECKABLE (A4)

Items for adjudication (manager notes, not rulings):

1. **A4 is UNCHECKABLE.** The failed r7 run never saved a cache. Its writeup (`r7_build_c2_false_positive.md`) records no built/None counts.
   - The one node-level fact it does record matches: (6.0, 2.8943) BUILT at 6.0086.
   - This run has 121 None / 3399 built. The per-node set is in `...built_none_set.json` on the Mac.
   - A full check needs the failed run's console log from the old build machine, if it still exists.
2. **G2:** max |zb residual| is 0.04981 km, which is under the 0.125 gate but over the r7 expectation of ≤0.035. 0.0498 sits just below the solver acceptance band solve_tol_km = 0.4 × 0.125 = 0.05 km. That is the early-exit mechanism ruled in §0.31.
3. **G5a, ocean mass residual:** one-signed as designed, from −0.00508 to −0.00037 (a 0.04–0.5% mass deficit). No ratified numeric bound exists.
4. **G3:** 0 `no_ocean_frozen_edge` Nones. All 121 Nones are C1-predicted reachability exclusions, so f_excluded = 0.034375 is entirely reachability.
5. **Merge metadata:** `build_mode` says "SLURM sharded (11 array tasks)", but the build ran as local processes. The label is hardcoded in the merge script.
6. **Cost:** 48.3 shard-hours, versus ~10 h serial for B's 2026-08-19 run. Per node this is about 5× slower, consistent with 11 concurrent processes on the Mac. It is not a correctness signal.

## Merge-script gate report (verbatim excerpt)

```
1_zb_placement_rejects_zero: 0 PASS
2_max_residual_under_zb_tol: 0.049808020610011994 km (threshold 0.125, expectation 0.035) PASS
3_reachability_none_set: n_ocean_none 121 / 3520, f_excluded 0.034375 PASS
4_moi_window_floor_remeasured: Cmeasured 0.335 ±0.08, n_nodes_measured 3399, max_abs_cmr2_deviation 0.0038613362093641856, margin 0.07613866379063582
5_frozen_invariants_not_fired: 39 built, 0 rejected, 0 if2_rejected PASS
6_f_excluded_recorded: 0.034375 PASS
total_shard_wall_h: 48.28101739784082
ACCEPTANCE GATES: ALL PASS
```

## enceladus_zbw_acceptance.py (verbatim)

```
[PASS       ] 0.30_zb_axis_88_nodes_end_45.0: n_zb=88, first=5.0, last=45.0
[PASS       ] r7_w_axis_40_nodes: n_w=40
[PASS       ] 0.30_ocean_nodes_3520: len(structures)=3520
[PASS       ] r7_frozen_nodes_39_axis: n_frozen=39
[PASS       ] r7_zb_tol_km_0.125: zb_tol_km=0.125
[PASS       ] config_declares_88x40: config n_total/n_w=(88, 40)
[PASS       ] A0_cache_table_equals_config_table: rows=40
[PASS       ] A0_band_is_0.4_zb_tol: reachability_band_km=0.05; predicate='excluded iff zb_target < zb_min(w) - reachability_band_km'
[PASS       ] A1_zero_C2_violations: 0
[PASS       ] A2_node_6.0_2.8943_built_zb_actual_6.0086: built=True, zb_actual=6.008573897790949
[PASS       ] A3_node_12.0_4.9239_None_unreachable: None=True, reason='build_exception', floor=12.0671, below-floor=0.0671 km (band 0.05)
[UNCHECKABLE] A4_built_none_set_identical_to_failed_r7_log: 121 None / 3399 built
[PASS       ] G1_zero_placement_rejects: 0
[PASS       ] G1b_placement_reject_reasons_zero: 0
[PASS       ] G2_max_residual_lt_0.125: max|residual|=0.04981 km; expectation (<=0.035) NOT met (informational)
[PASS       ] G3_none_set_equals_C1_prediction: 0 nodes outside prediction; 121 total None, 0 no_ocean_frozen_edge
[PASS       ] G4_moi_floor_remeasured: recorded == recomputed 0.0038613362093641856, margin 0.0761, n_measured 3399
[PASS       ] G5a_ocean_mass_invariant_every_node: n=3399, min=-0.005079502934320468, max=-0.0003722761790641016
[PASS       ] G5b_frozen_mass_invariant_I-F1_every_node: n=39, max|res|=3.33e-16 (<= 1e-06)
[PASS       ] G6a_frozen_no_rejects_I-F1..I-F4: all rejects 0, n_built=39
[PASS       ] G6b_frozen_nodes_recheck: max closure 9.09e-13 kg/m3 (<=12.0), max|zb resid| 0.0978 km (<0.25), liquid layers all 0
[PASS       ] G6c_frozen_I-F6_non_conditioning: 4/39 inside 0.335±0.001; span 0.32400-0.34392
[PASS       ] G7_f_excluded_recorded: 121/3520 = 0.034375; C3 shift +log(1/(1-f)) = 0.034980
[PASS       ] S_shards_complete_and_match_merged: 11 shards, mask equal, no quarantine
SUMMARY {'PASS': 23, 'FAIL': 0, 'UNCHECKABLE': 1}
```

## Scientific-reviewer adjudication (2026-10-21): ACCEPT WITH CONDITIONS

- **A4: WAIVED (deductive equivalence), not PASS.**
  - Recomputing both predicates over the 88×40 axes, the old predicate excludes 122 nodes and the new one 121. They differ only at (6.0, 2.8943), whose floor is 6.0088.
  - C2 is a post-build check (cache_builder.py:2144-2153), so the predicate does not change which nodes build. Between b7beb6b4 and 634ea0cd the only non-text change is docstring-only (LayerPropagators.py).
  - The failed run's single violation was that node. Its built set therefore equals the new non-excluded set: 3399 built / 121 None, and the predicted set difference is empty.
  - The nearest node to the band edge has a 0.017 km margin.
  - If the old console log turns up, diff it against this set and append the result here.
- **G2: benign.** solve_tol_km = 0.4 × 0.125 = 0.05 (cache_builder.py:1996) is passed as zbTol_km (:2053). IceShellResidual stops on ANY iterate within tolerance (LayerPropagators.py:667, caught at :733), so residuals are bounded just under 0.05 km.
- **C4:** the r7 expectation "≤0.035 km" is SUPERSEDED by the <0.05 km bound from the solver's early exit.
- **G3: consistent.** The frozen axis starts at 46.5 km, so no edge Nones are expected.
- **G5a: not blocking.** It lies inside the deficit band of already-ratified caches (Titan to −0.98%, Europa −0.28%; mass_gate_blast_radius_audit.json). However:
  - σ/|obs| is 0.68% for C20 and 0.97% for C22, so a 1:1 propagation of a −0.51% deficit into C20 would be ≈0.75σ.
  - The frozen branch closes mass exactly, so the deficit is a branch asymmetry that feeds log B.

Conditions:
- **C1:** every one of the 121 `build_exception` Nones must be the bracket / same-sign "target outside achievable zb" failure (no I/O, cache-path or EOS errors), checked from the shard logs with `plans/cluster/enceladus_zbw_none_reasons.sh`. A3 depends on this too.
- **C2:** merge report provenance. `build_mode` was hardcoded to SLURM; the merge script now records the real mode plus build_commit, merge_commit and merge_host. Re-run the merge to regenerate `r7_production_build_report.json`.
- **C3, before ANY evidence or observable is quoted:**
  - Rebuild the ~5 largest-deficit nodes with a finer silicate radial step.
  - The C20/C22/libration shift must be ≤0.1 σ_obs; otherwise bound or correct the deficit.
  - Map the deficit over (zb, w) to check for jagged noise.
  - The study script is `plans/scripts/enceladus_zbw_mass_deficit_study.py` (in preparation).
- **C4:** applied above.

Gating:
- The cache may move into `PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/` once C1 and C2 are met.
- The two-run reference protocol (§0.30 item 5) may START, but no branch evidence may be quoted until C3 passes.
