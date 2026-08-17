# r7 — RATIFIED WITH CONDITIONS (reviewer, 2026-08-16). BUILD AUTHORIZED after C1-C3.

E1(a)/(b)/(c) + E2 all RATIFIED (independent 12-node production build:
12/12, 0 rejects, max residual 27% of budget; E1(c) equality to 5e-12 km
on real Planet objects + convention-consistency + frozen-branch-parity
reasons; bracket reproduces vs a cache-free continuum to 3e-7 km).
Suite 319+1xfail+2 pre-existing confirmed.

PRE-BUILD BLOCKING:
C1: map the ocean reachability restriction over ALL 40 w nodes (procedure:
 one build_single_structure at zb=1 km per w with zbTol_km set; read
 zb_min/max off the bracket-end evaluations; ~15 min) + config block
 metadata.structure_cache_spec.ocean_reachability_restriction with the
 40-point table, mechanism (GetTfreeze 0.05 K step), ±6.3 km boundary
 uncertainty, excluded fraction 8.9-24% band, "model-reachability
 restriction NOT a physical support edge" in those words, and the
 no-posterior-mass cross-reference (region below ~18 km, >=8 sigma
 rejected; non-monotone in w, peak near 5 ppt: zb_min 5.64/9.18/11.55/
 7.43/8.42/9.49/2.84 km at w=0.1/1/5/10/20/50/100).
C2: builder records the declared restriction and FAILS on any None
 ocean node outside it (two-sided: must reject inside, build outside).
C3: evidence normalization +log(1/(1-f_excluded)) from the DELIVERED
 cache's None mask recorded in branch_model.evidence_protocol as a
 validity condition BEFORE any log B; folds into B-A2 (trigger fired:
 f>5%); report MoI-window and reachability exclusions as separate lines.
PRE-FREEZE NON-BLOCKING: C4 byte-identity wording/test-name corrections
 (legacy path returns best-iterate too — campaigns safe: ICEIh_THICKNESS
 unused by any ratified cache; inductogram/GUI improve); C5 undercount
 figure is zb-dependent 0.26-0.86 km (~zb/52), not 0.44-0.50; C6 state
 the E2 bracket check used synthetic nodes — B5 NOT discharged.

MACHINE B BUILD INSTRUCTIONS (after C1-C3; full text in this file's
source adjudication): build_zbw_grid_cache with PPEnceladus template;
zb segmented axis 5-20@1.0 / 20-22@0.5 / 22-30@0.25 / 30-42@0.5 /
42-45@0.25 = 87 nodes; w = 10**linspace(-1,2,40); ocean_overrides
{comp: Seawater, deltaT: 0.002}; zb_tol_km=0.125 EXPLICIT; config=<the
candidate> and NO bulk_overrides argument (provenance chain: assert
ocean_moi_window source == config metadata, Cuncertainty 0.08/0.08);
frozen axis arange(46.5, 65.8001, 0.5)=39 nodes, frozen_Cuncertainty
0.015, mass_tol 1e-6, rho_closure 12.0, zb_tol 0.25,
moi_nonconditioning_window None (I-F6 tests template 0.335±0.001),
max_iter 10; extrap_ocean False; tb_placeholder 272.0. Schema
v3.2-zbw-joint. ACCEPTANCE: 0 placement rejects, max residual <0.125
(expect <=0.035); None set == C1 prediction within one cell (any outside
= halt); re-measured MoI floor at production scale; mass invariant on
every node; frozen I-F1..I-F4 + I-F6 must NOT fire; f_excluded recorded.
Budget ~20 s/node: ocean ~19.5 h serial (parallelize by segment), frozen
~35 min. THEN: two nested-sampling runs (ocean 5-param, frozen 3-param),
n_eff=8000 n_active=4096, 3 seeds each, no branch coordinate; log B post
hoc at 0.5/0.5 with the C3 normalization; one-sided if underflow; raise
n_eff rather than quote sign if |log B| < 3x numerical error. B-A1
(rock-model delta, <=0.05 sigma on all rows voids asymmetry 1) and B-A2
(realized cmr2 + BOTH exclusions separately) in the run design.
