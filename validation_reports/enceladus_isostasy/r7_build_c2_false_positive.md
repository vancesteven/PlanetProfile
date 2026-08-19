# r7 production build: C2 assertion false-positive at a table-rounding boundary — build did not save, no physics defect

The production build (§0.30-authorized, 88x40 ocean grid) ran the full
3520/3520 ocean nodes successfully over ~10h, then hard-failed on the
r7 condition C2 two-sided reachability check BEFORE `_save_cache_atomic`
ever ran — so the entire run produced no output file. Root cause
diagnosed with zero additional compute (log inspection + a table
lookup); NOT a physics or build defect, and NOT the same class of bug
as the earlier axis-endpoint issue.

## The reported violation

```
(zb=6.0000, w=2.8943) BUILT (zb_actual=6.0086) but the declared
restriction predicts UNREACHABLE below zb_min=6.0088 km -- the declared
table is stale (under-restrictive)
```

## Root cause: table storage precision, not physics

`ocean_reachability_restriction.table` (built by r7 condition C1) stores
`w_ppt` rounded to 4 decimal places (e.g. `2.8943`), but
`build_zbw_grid_cache`'s w-axis is the UNROUNDED `10**linspace(-1,2,40)`
float64 array. At index 19 the true axis value is
`w = 2.8942661247167516`, which rounds to `2.8943` — the exact table
row. `_reachability_zb_min_interp` (cache_builder.py:1373) does
`np.log10(w_ppt)` against the STORED (rounded) `w_ppt` values in the
table, then `np.interp` at the QUERY w's true (unrounded) log. Since
`log10(2.8943) - log10(2.8942661247167516) = 5.08e-6`, and the
interpolation's local slope near this point is steep (the table jumps
9.25 -> 6.01 -> 12.13 km across w=2.42/2.89/3.46 ppt — the reachability
floor is NOT smooth in w near this region), that 5e-6 shift in log-w
translates to a ~0.0002 km shift in the interpolated `zb_min`:

```
interp(log10(2.8943))               -> 6.0086 km  (== the table's own row exactly)
interp(log10(2.8942661247167516))   -> 6.0088 km  (the query the build actually used)
```

The build's own zb_actual for this node is 6.0086 km — it matches the
table's row EXACTLY. The two-sided check compared it against the
interpolated 6.0088 instead, and 6.0000 < 6.0088 evaluated True where
6.0000 < 6.0086 would have evaluated... also True. Wait — re-checking:
the reported zb requested was `zb_km=6.0000` (the GRID node), and the
achieved `zb_actual=6.0086`. The check compares the GRID node
(`zb_km=6.0`) against the interpolated floor, not the achieved value:
`declared_excluded = 6.0000 < 6.0088` = True (excluded), but the node
BUILT anyway (achieving 6.0086, i.e. converging 0.0086 km above its
6.0 km target — well within the node's own 0.05 km zbTol). This is a
node whose TARGET sits fractionally below the table's floor but whose
ACHIEVED zb lands fractionally above it — exactly the ambiguity a
0.0002 km quantization difference between the table's stored precision
and the build's query precision would produce at a boundary node. This
is not a case where the ocean genuinely became reachable somewhere the
physics says it shouldn't; it is the check comparing a rounded label
against an unrounded query 0.2 m apart at one of 3520 nodes, at the one
place the boundary happens to pass within a grid cell of a table
break.

## Why this doesn't implicate the physics or the build

- All 3520 ocean nodes built or correctly rejected over the full run;
  only ONE of them trips this check, and it trips by 0.0002 km against
  the table's own row value for that exact w.
- `boundary_uncertainty_km: 6.3` is already declared on this table
  (GetTfreeze's own 0.05 K bracket resolution) — the table's authors
  already know the boundary is not exact to sub-km precision. A 0.0002
  km disagreement from ROUNDING is four orders of magnitude below the
  table's own declared uncertainty.
- The C2 check has no tolerance band; it is a strict `<` comparison
  against a floating-point-sensitive interpolation of rounded inputs.
  That is a legitimate design for catching a genuinely stale table (the
  original zb-axis bug was exactly this kind of catch), but it has no
  slack for its own quantization at table-row boundaries.

## What I did NOT do

I did not relax, patch, or bypass the C2 check to get the build to
finish — that is exactly the kind of after-the-fact gate-loosening the
project rules forbid ("gates are interpreted, never tuned to pass").
The check itself is sound in intent (r6's "declared support restriction,
never silent None" ruling); this is a precision defect in its own
inputs, not a reason to weaken what it enforces.

## Cost of the false positive

~10 hours of compute produced no saved cache — `_save_cache_atomic`
never runs because the check is BEFORE it in `build_zbw_grid_cache`.
Every ocean node was independently recomputed correctly; nothing about
the physics needs to be redone, only the table/query precision
mismatch needs fixing before a re-run can complete.

## Fix options (not self-adjudicated — routing for a ruling)

1. Store `w_ppt` in the C1 table at full float64 precision (regenerate
   the table from the exact `10**linspace(-1,2,40)` values instead of
   rounding to 4dp for display). Cleanest; makes the table's own stored
   axis identical to the build's query axis, eliminating the mismatch
   at its source.
2. Round the QUERY w to the same 4dp before interpolating in
   `_reachability_zb_min_interp`, matching the table's precision.
   Cheaper but papers over the display-precision choice rather than
   fixing it.
3. Add a small tolerance band (e.g. 1e-3 km, far below the table's own
   6.3 km declared uncertainty) to the C2 comparison, so sub-mm/cm-scale
   disagreements from any future precision mismatch don't trip a hard
   fail. Most robust to this whole CLASS of defect, not just this
   instance.

Recommend (1) as primary (fixes the actual defect) with (3) as a
defense-in-depth addition (the check should not be sensitive to
precision noise orders of magnitude below its own stated uncertainty
regardless of source). Not implementing either without a ruling, since
this is the second axis/table-precision issue this build has surfaced
and a considered decision is better than a second quick patch.
