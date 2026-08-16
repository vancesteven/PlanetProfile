# r6 — NOT RATIFIED (reviewer, 2026-08-16). BUILD AUTHORIZATION WITHHELD.

D1-D8 repairs all reproduce adversarially (bracket/floor/deltas exact).
TWO NEW STRUCTURAL BLOCKERS, both invisible to bypass-harness testing:

E1 (BUILD): zb-placement invariant fails GRID-WIDE — 8/8 real builds
rejected by 5-27x (residuals up to -3.4 km at the posterior heart).
Root cause: GetIceShellTFreeze root-finds Tb at xtol=TfreezeRes_K=0.05 K
while the whole 5-45 km range spans ~0.27 K (0.05 K ~ 7 km of zb), and
LayerPropagators best_planet (585-587) retains the LAST iterate, not the
best. Every campaign test used zb_tol_km=5.0 (40x production); the D1
probe bypassed the invariant via build_single_structure. Production
build would abort with an empty cache. Also contaminates D1's 8-node
numerical evidence (structural nestedness argument STANDS; 0.08 stands;
re-measure the floor on admissible nodes).
E2 (TRAINING): no ocean-branch (zb,w) cache lookup exists —
_struct_for_hydrosphere falls through and returns the WHOLE CACHE DICT;
every ocean sample rejects to -inf. The frozen analogue was written
(A8); the ocean one never was. D2's check ran on a Tb-keyed smoke cache.

Rulings: D1 structural argument + 0.08 RATIFIED (numerical half to
re-measure); branch-scope transfer RATIFIED as a consistency convention
(wording condition); D8 disposition RATIFIED but justification REFUTED
(exact projection IMPROVES the B13 gate, L 1.91->1.19, same solution —
favors option (i); fix the budget phrasing incl. the one-sided-bias-as-
symmetric-sigma clause); corner-nodes item MOOT (failure is grid-wide;
re-pose post-E1; rule then = declared support restriction, never silent
None). Park erratum gap under-recorded; PDF still unfiled (user).

r7 path: 1) E1 — tighten Tb root-find + fix best_planet to retain the
minimum-residual iterate; verify at PRODUCTION tolerance on >=1 node per
segment; 2) E2 — ocean (zb,w) lookup mirroring _frozen_struct_for_theta;
3) re-measure D1 floor on admissible nodes; 4) D8/Sigma_model phrasing +
erratum gap spec; 5) re-run B1; r7.

STANDING RULE (manager-adopted): ALL verification runs through the
production path at production tolerances — no build_single_structure
bypasses, no relaxed zb_tol_km, no smoke-cache stand-ins. The r5 "through
the production dispatch" instruction generalizes to cache + dataset paths.
