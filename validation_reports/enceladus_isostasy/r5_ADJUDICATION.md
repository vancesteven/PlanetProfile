# r5 ratification pass — NOT RATIFIED (reviewer, 2026-08-16)

Suite 86+1xfail verified; bijection/support/CMR2-non-conditioning/
Sigma_model/Park values/H&M constants/A8 dispatch all reproduce exactly.
Grid segmentation robust to the pending repairs. Provenance gap: Park
erratum PDF not in papers/ (repo copy is the Jan-2024 original; 0.33
sigma difference on libration) — file it.

BUILD BLOCKERS:
- D1: ocean-branch MoI window unspecified/unenforceable — bulk_overrides
  is not a config field; builder defaults to template Cuncertainty
  0.001; the asserted ±0.08 has no recorded derivation anywhere (Titan
  used 0.06 via explicit overrides). = MODERATE-4 ocean half, never
  closed.
- D2: generate_sbi_dataset cannot produce the observable vector under
  isostatic_hm2019 (gravity_pair None -> C20/C22 fall through; NO C30
  arm; drop_nonfinite rejects everything -> empty dataset). C30 was
  wired for the likelihood only.

LIKELIHOOD DEFECTS (all on the sole thickness channel; repair TOGETHER,
verify through the production dispatch, not librations() directly):
- D3: the B2'-ruled Delta_rho treatment is NOT wired at either
  production call site (mcmc_runner 2428/2486 run the historical
  hydrostatic convention; 1.3-2.6 km on the headline deliverable;
  dlibration/dC2 identically zero as shipped).
- D4: libration_model_correction (1.008, +0.245 sigma) never applied.
- D5: libration_sys_frac inert (sampled, never in chi2) — dead param by
  the campaign's own CRITICAL-1 criterion; the ±0.4% residual it
  marginalizes is unmarginalized.
- D6: rho_ice_kgm3 never reaches the ocean-branch libration
  (rho_shell_override unreachable from the likelihood) — a THIRD branch
  asymmetry, in the likelihood; REPAIR, not bound (B-A3).
- D7: rho_rock_kgm3 missing from PARAMETER_REGISTRY.
- D8: first-order H22 term (0.65 sigma one-sided on C22) unbudgeted
  (Sigma_model C22 is 5.6x smaller).

ASYMMETRY RULING (binding): BOUND AND CAVEAT. Key finding: the porous
profile never reaches the conditioned observables (both branches reduce
to volume-weighted uniform interiors) — asymmetry (1) collapses into
structure selection. Equalizing rejected both directions. Required in
B's run design: B-A1 rock-model delta measurement (<=0.05 sigma on all
rows -> asymmetry VOID); B-A2 realized-cmr2 distribution + excluded
prior-volume fraction (>5% -> odds reported over the window-restricted
prior); B-A3 = D6 repair. Verdict robust (23-45 sigma gravity, 19-20
sigma libration); caveats govern the published number only.

FREEZE-GATE CORRECTION (binding): Tajeddine 0.08 sigma is ALREADY in
Sigma_model (libration_deg 0.00025) — must not also enter the
libration_sys_frac width (double-count). det(A_mat) margin 0.9 km
recorded, still unimplemented. Freeze gate moot until D5.

TRIAGE: closed = B4/schema/Sigma_model/C30-likelihood/registry(4 of 5).
Build-blocked = D1, D2. Training-blocked = D3-D6. Publication-blocked =
B1 (re-run after D3), B5, B12, B8/B9 display, D8. Plus listed config
documentation inconsistencies (incl. a quote of the forbidden rescinded
25.80 km number in the segment rationale).

RE-SUBMISSION PATH (ordered): 1) settle D1 (cache build is provably
independent of D3-D6 — cache_builder has zero libration/isostasy refs —
so the BUILD can start the moment D1 closes, parallel to the repairs);
2) repair D3-D6 together + numerical verification reproducing the ruled
bracket through the production dispatch; 3) wire D2 before dataset gen;
4) re-run B1 on the corrected model (expect ~2-2.5 km shift), disposition
D7/D8/doc items, re-submit for r6.
